<script type="text/javascript">
////////////////////////////////////////////////////////////////////////////////
// Data     																///
//////////////////////////////////////////////////////////////////////////////
// according to OpenSMILES (http://opensmiles.org/opensmiles.html)
const aromatics_map = new Map();
const aromatics_set = new Set(['b', 'c', 'n', 'o', 'p', 's', 'as', 'se', 'te']);
const organic_set = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'b', 'c', 'n', 'o', 'p', 's']);
const bonds_set = new Set(['-', '=', '#', '$', ':', '\\', '/']);
//Allowed characters based on the analysis of SMILES from ChEMBL v35, excluding "." - multicomponent structures are not allowed.
const smiles_chars = new Set(["C", "c", "1", "(", "-", "n", "2", "=", "O", ")", "[", "H", "]", "l", "#", "N", "3", "B", "r", "S", "+", "/", "I", "\\", "4", "s", "5", "F", "o", "@", "6", "7", "8", "9", "%", "0", ".", "P", "a", "i", "L", "K", "e", "t", "T", "Z", "M", "g", "A", "p", "R", "b", "X", "G"]);
//SEE: https://en.wikipedia.org/wiki/List_of_chemical_elements
const atom_symbols = new Set(["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "b", "c", "n", "o", "s", "p", "te", "as", "se"]);
const bonds_map = new Map();
//Add valences to aromatics
aromatics_map.set("b", new Set([3]));
aromatics_map.set("c", new Set([4]));
aromatics_map.set("n", new Set([3,5]));
aromatics_map.set("o", new Set([2]));
aromatics_map.set("p", new Set([3,5]));
aromatics_map.set("s", new Set([2,4,6]));
//Add `valences` to bonds
bonds_map.set("-", [1]);
bonds_map.set("=", [2]);
bonds_map.set("#", [3]);
bonds_map.set("$", [4]);
bonds_map.set(":", [1]);
bonds_map.set("\\", [1]);
bonds_map.set("/", [1]);



////////////////////////////////////////////////////////////////////////////////
// Functions to deal with the atoms and bonds								///
//////////////////////////////////////////////////////////////////////////////


// 1	Remove stereo
function clear_stereo (smiles) {
	smiles = smiles.replaceAll("@TB20", "");
	smiles = smiles.replaceAll("@OH30", "");
	smiles = smiles.replaceAll(/@TB1[0-9]/g, "");
	smiles = smiles.replaceAll(/@OH1[0-9]/g, "");
	smiles = smiles.replaceAll(/@OH2[0-9]/g, "");
	smiles = smiles.replaceAll(/@TB[1-9]/g, "");
	smiles = smiles.replaceAll(/@OH[1-9]/g, "");
	smiles = smiles.replaceAll("@TH1", "");
	smiles = smiles.replaceAll("@AL1", "");
	smiles = smiles.replaceAll("@AL2", "");
	smiles = smiles.replaceAll("@SP1", "");
	smiles = smiles.replaceAll("@SP2", "");
	smiles = smiles.replaceAll("@SP3", "");
	smiles = smiles.replaceAll("@@", "");
	smiles = smiles.replaceAll("@", "");
	return(smiles);
}

// 2 Decide on the isotope of this atom
function get_isotope (atom_str, position) {
	// Position is in the SMILES substring corresponding to this atom
	let isotope = "";
	let is_isotope = 0;
	if ((/[0-9]/).test(atom_str[position])) {
		is_isotope = 1;
		isotope = isotope + atom_str[position];
		position++
		//Get the whole isotope label
		while ((/[0-9]/).test(atom_str[position])) {
			isotope = isotope + atom_str[position];
			position++
		}
	}
	return {"position": position, "is_isotope": is_isotope, "isotope": isotope};
}

// 3 Get the SMILES substring corresponding to this bracketed atom 
function get_symbol_bracket (atom_str, position) {
	let atom_symbol = "";
	while ( (/[a-zA-Z]/).test(atom_str[position]) ) {
		atom_symbol = atom_symbol + atom_str[position];
		position++
		if (position >= atom_str.length) {
			return {"msg": 'Unable to parse part of the input: ' + '\r\n' + atom_str + '...\r\n' + 'Please, make corrections or consider alternative input'};		
		}
	}
	return {"position": position, "symbol": atom_symbol};
}

// 4 Get H-status and H-count for this atom
function get_h (atom_str, position) {
	// Position is the position of the first H belonging to this atom in this atom SMILES substring
	let explH = 0;
	let explH_count = 0;
	if (atom_str[position] === "H") {
		position++;
		explH = 1;
		if ( (/[0-9]/).test(atom_str[position]) ) {
			explH_count = atom_str[position];
			position++
		} else {
			explH_count = explH;
		}
	}
	return {"position": position, "explH": explH, "explH_count": explH_count};
}

// 5 Get the charge of this atom
function get_charge (atom_str, position) {
	// Position is in the SMILES substring corresponding to this atom
	let is_charged = 0;
	let charge_qual = "";
	let charge_quant = 0;
	if ( (/\-|\+/).test(atom_str[position]) ) {
		is_charged = 1;
		charge_qual = atom_str[position];
		charge_quant = 1;
		position++;
		//Check for quantifier
		if (atom_str[position] === charge_qual) {
			charge_quant = 2;
			position++;
		} else {
			//Digital quantifier
			if ((/[0-9]/).test(atom_str[position])) {
				if ((/[0-9]/).test(atom_str[position + 1])) {
					charge_quant = atom_str[position] + atom_str[position + 1];
					position++;
					position++;
				} else {
					charge_quant = atom_str[position];
					position++
				}
			}
		}
	}
	return {"position": position, "is_charged": is_charged, "charge_qual": charge_qual, "charge_quant": charge_quant};
}

// 6 Get the bonds on the left for this atom
function pre_bond_get (position, atom_start, smiles, atom_bonds, bonds_set) {
	// Position is the current position in the whole SMILES string
	let found = 0;
	let bond = "";
	if (bonds_set.has(smiles[atom_start-1])) {
		found ++
		atom_start = atom_start-1;
		bond = smiles[atom_start];
		//Problem, report:
		if (atom_start === 0) {
			return {"msg": 'Unable to parse input starting with bond: ' + '\r\n' + smiles.substring(0, position) + '...\r\n' + 'Please, make corrections or consider alternative input'};
		}
		return {"found": found, "atom_start": atom_start, "bond": bond};
	}
}

// 7 Get the connected atom on the left if the previous atom in the SMILES string is in branch
function prev_atom_get_main (position, atom_start, smiles) {
	// Position is the current position in the whole SMILES string
	let connected_position = -1;
	if (smiles[atom_start-1] === ")") {
		let open_branch = 1;
		for (let k = atom_start-2; k >= 0; k--) {
			//Count '(' / ')'
			if (smiles[k] === ')') {
				open_branch++;
			}
			if (smiles[k] === '(') {
				open_branch--;
			}
			k--
			if (smiles[k] === ')') {
				open_branch++;
			}
			if (smiles[k] === '(') {
				open_branch--;
			}
			if (open_branch === 0) {
				connected_position = k-1;
				break;
			}
			if (k === 0) {
				return {"msg"	:'Unable to find the opening bracket "(" not in the first position in: ' + '\r\n' + smiles.substring(0, position) + '\r\n' + 'Please, make corrections or consider alternative input'};
			}
		}
	return {"position": position, "atom_start": atom_start, "connected_position": connected_position};	
	} else {
	return "not_the_case";
	}
}

// 8 Get the connected atom on the left if this atom is the start of the branch
function prev_atom_get_frombranch (position, atom_start, smiles) {
	// Position is the current position in the whole SMILES string
	let connected_position = -1;
	if ( smiles[atom_start-1] === "(" ) {
		atom_start--
		let pointer = atom_start-1;
		for (let k = pointer; k >= 0; k--) {
			if (k === 0) {
				return {"msg": "Unable to parse input starting with '(': " + "\r\n" + smiles.substring(0, position) + "...\r\n" + "Please, make corrections or consider alternative input"}
			}
			if (smiles[k] === "(") {
				pointer = k;
				continue;
			} else {
				pointer = k;
				break;
			}
		}
		if (smiles[pointer] != ")") {
			connected_position = pointer;
			return {"position": position, "atom_start": atom_start, "connected_position": connected_position};
		} else {
			let open_branch = 1;
			for (let k = pointer-1; k >= 0; k--) {
				//Count '(' / ')'
				if (smiles[k] === ')') {
					open_branch++;
				}
				if (smiles[k] === '(') {
					open_branch--;
				}
				k--
				if (smiles[k] === ')') {
					open_branch++;
				}
				if (smiles[k] === '(') {
					open_branch--;
				}
				if (open_branch === 0) {
					connected_position = k-1;
					break;
				}
				if (k === 0) {
					return {"msg":"Unable to find the opening bracket '(' not in the first position in: " + "\r\n" + smiles.substring(0, position) + "\r\n" + "Please, make corrections or consider alternative input"};
				}
			}
			return {"position": position, "atom_start": atom_start, "connected_position": connected_position};
		}
	} else {
		return "not_the_case";
	}
}

// 9 Get the connected atom on the left if the previous symbol is digit
function prev_atom_get_digit (atom_start, smiles) {
	let connected_position;
	if ( (/[0-9]/).test(smiles[atom_start-1]) ) {
		connected_position = atom_start-1;
		return {"atom_start": atom_start, "connected_position": connected_position};
	} else {
		return "not_the_case";
	}
}

// 10 Get the connected atom on the left if the previous symbol is exactly atom
function prev_atom_get_atom (atom_start, smiles) {
	let connected_position;
	if ( organic_set.has(smiles[atom_start-1]) | smiles[atom_start-1] == "]" ) {
		connected_position = atom_start-1;
		return {"atom_start": atom_start, "connected_position": connected_position};
	} else {
		return "not_the_case";
	}
}

// 11 Get ring info from the right side of the SMILES string
function next_get_rings (position, atom_end, smiles) {
	// Position is the current position in the whole SMILES string
	position++
	let ring_str = "";
	let ring_positions = new Set();
	let ring_start;
	if ( (/[0-9]|\%/).test(smiles[position]) ) {
		ring_start = position;
		ring_str = ring_str + smiles[position];
		ring_positions.add(position);
		atom_end = position;
		position++
		while ( (/[0-9]|'%/).test(smiles[position]) ) {
			ring_str = ring_str + smiles[position];
			ring_positions.add(position);
			atom_end = position;
			position++
		}
	}
	return {"position": position, "atom_end": atom_end, "ring_str": ring_str, "ring_positions": ring_positions, "ring_start": ring_start};
}

// 12 Parse the ring data into the distinct rings
function parse_ring_str (ring_start, ring_str) {
	let rings = [];
	for (let i = 0; i < ring_str.length; i++) {
		if (ring_str[i] === '%') {
			let this_ring = {};
			let this_ring_str = ring_str[i] + ring_str[i+1] + ring_str[i+2];
			let this_ring_pos = new Set([ring_start + i, ring_start + i + 1, ring_start + i + 2]);
			this_ring.str = this_ring_str;
			this_ring.pos = this_ring_pos;
			rings.push(this_ring);
			i = i + 3;
		}
		if ( (/[0-9]/).test(ring_str[i]) ) {
			let this_ring = {};
			let this_ring_str = ring_str[i];
			let this_ring_pos = new Set([ring_start + i]);
			this_ring.str = this_ring_str;
			this_ring.pos = this_ring_pos;
			rings.push(this_ring);
		}
	}
	return rings;
}

// 13 Get the ending positon of this ring
function end_ring (atom_end, this_ring, smiles) {
	let ring_symbol = this_ring.str;
	let ring_pattern = new RegExp(ring_symbol + "(?![^\\[*]})", "g");
	let ring_end = smiles.substring(atom_end+1, smiles.length).search(ring_pattern);
	if (ring_end > -1) {
		return {"end_start": ring_end + atom_end + 1, "end_end": ring_end + atom_end + 1 + ring_symbol.length -1};
	} else {
		return -1;
	}
}


////////////////////////////////////////////////////////////////////////////////
// Functions to do the basic check-up										///
//////////////////////////////////////////////////////////////////////////////
// What could be wrong with the SMILES string provided by the user?
// - Length is to big
// - Line breaks, multiline input
// - Forbidden characters

// 1 Validate input string
function validate_smiles_input (smiles_input) {
	let smiles;
	let cr_eol;
	let lf_eol;
	let space_dlm;
	let tab_dlm;
	// Trim the input first
	smiles = smiles_input.trim();
	// Locate the first EOL, which will mean that input is multilined
	// If so, -> delete * starting from EOL and trim again just in case
	cr_eol = smiles.indexOf('\r');
	if (cr_eol > -1) {
		smiles.slice(0, cr_eol).trim();
	}
	lf_eol = smiles.indexOf('\n');
	if (smiles.indexOf(lf_eol) > -1) {
		smiles.slice(0, lf_eol).trim();
	}
	// Locate the first space or tabulation, which will mean that string contains several chemical structures
	// Process only the first one
	space_dlm = smiles.indexOf(' ');
	if (space_dlm > -1) {
		smiles.slice(0, space_dlm).trim();
	}
	tab_dlm = smiles.indexOf('	');
	if (tab_dlm > -1) {
		smiles.slice(0, tab_dlm).trim();
	}
	// Check the length
	if (smiles.length > 500) {
		return("smth_wrong");
		alert("Input string is too long");
	}
	if (smiles.length < 3) {
		alert("Input string is too short");
		return("smth_wrong");
	}
	// Return the SMILES string
	return(smiles);
}

// 2 Check input smiles for components
function select_component (input_smiles) {
	let smiles;
	let point_dlm;
	point_dlm = input_smiles.indexOf('.');
	if (point_dlm > -1) {
		let smiles_array = input_smiles.split('.');
		smiles = smiles_array.reduce( function (a, b) { return a.length > b.length ? a : b; } ).trim();
	}
	// Check the length
	if (smiles.length > 500) {
		return("smth_wrong");
		alert("Input string is too long");
	}
	if (smiles.length < 3) {
		alert("Input string is too short");
		return("smth_wrong");
	}
	// Return the SMILES string
	return(smiles);
}

////////////////////////////////////////////////////////////////////////////////
// Functions to process SMILES												///
//////////////////////////////////////////////////////////////////////////////

// Remove flags on stereochemistry, and check length again
function destereo_smiles (smiles_input) {
	const smiles = clear_stereo(smiles_input);
	if (smiles.length < 3) {
		alert("Input SMILES string is too short");
		return("smth_wrong");
	}
	return(smiles);
}

// Parse SMILES,
// #red this function is still immature and should checked latter: some unused vars could be present and some exceptions could be unsecured, some variables have questionable names, some loging is no longer needed
function parse_smiles (smiles, aromatics_map, aromatics_set, organic_set, bonds_set, smiles_chars, atom_symbols, bonds_map) {
	//create data structures to store the structure
	mol_structure_raw = [];
	//create set to store the positions of the already connected rings
	let closed_rings = new Set();
	let atom_id = -1;
	//Traverse the SMILES from Left to Right
	smiles_traversal: for (let i = 0; i < smiles.length; i++) {
		//Create variable to distinguish between the bracketed atoms and atoms from organic subset
		let bracketed = 0;
		let organic = 0;
		//Declare the variables to store the characteristics of the atom
		let atom;
		let atom_start = -1;
		let atom_end = -1;
		// Basic descriptors
		let explH;
		let explH_count;
		let is_charged;
		let is_aromatic;
		let charge_qual;
		let charge_quant;
		let is_isotope;
		let isotope;
		let atom_symbol;
		let atom_bonds;
		// Store connections
		let connected_positions;
		//Check if current position belongs to the bracketed atom and characterize it
		if ( smiles[i] === '[') {
			console.log("Atom:  " + i);
			bracketed = 1;
			//Initialize the variables to store the characteristics of the atom
			atom = {};
			// Basic descriptors
			explH = 0;
			explH_count = 0;
			is_charged = 0;
			is_aromatic = 0;
			charge_qual = "";
			charge_quant = 0;
			is_isotope = 0;
			isotope = "";
			atom_symbol = "";
			atom_bonds = [];
			// Store connections
			let connected_positions = [];
			//positions
			atom_start = i;
			atom_end = i;
			//New atom
			atom_id++;
			atom.atom_id = atom_id;
			atom.bonds = [];
			
			//Get the whole atom substr
			let pos_end = i;
			close_sqbracket: for (let k = i+1; k < smiles.length; k++) {
				if (smiles[k] === ']') {
					pos_end++
					break close_sqbracket;
				}
				if (k === smiles.length-1 && smiles[k] != ']') {
					alert("] is missing");
					break smiles_traversal;
				}
				pos_end++
			}
			//Get the substring corresponding to this section
			let atom_str = smiles.substring(i, pos_end + 1);
			//Traversing the atom str leaving the opening bracket behind
			let atom_i = 1;
			//Get the isotope
			let isotope_prop = get_isotope(atom_str, atom_i);
			if (isotope_prop.is_isotope === 1) {
				//Add the data
				atom.is_isotope = 1;
				atom.isotope = isotope_prop.isotope;
				//Increase counter
				atom_i = isotope_prop['position'];
			} else {
				atom.is_isotope = 0;
				atom.isotope = "";
			}
			//Get atom symbol
			atom_symbol = get_symbol_bracket(atom_str, atom_i);
			if (typeof atom_symbol['msg'] === "undefined") {
				atom_i = atom_symbol['position'];
				atom_symbol = atom_symbol['symbol'];
			} else {
				alert(atom_symbol['msg']);
				break smiles_traversal;
			}
			//Check atom symbol and delete H if needed
			if(!atom_symbols.has(atom_symbol)) {
				if (atom_symbol.slice(-1) === 'H') {
					//Atom has explicit hydrogens, delete 'H' and decrease counter
					atom_symbol = atom_symbol.slice(0, atom_symbol.length-1);
					atom_i = atom_i-1;
					//Check symbol again
					if(!atom_symbols.has(atom_symbol)) {
						alert('Unknown atom symbol: ' + '\r\n' + atom_symbol + '\r\n' + 'Please, make corrections or consider alternative input');
						break smiles_traversal;
					}
				} else {
					alert('Unknown atom symbol: ' + '\r\n' + atom_symbol + '\r\n' + 'Please, make corrections or consider alternative input');
					break smiles_traversal;
				}
			}
			//Add to the description of atom
			atom.atom_symbol = atom_symbol;
			//Check if aromatic
			if (aromatics_set.has(atom_symbol)) {
				is_aromatic = 1;
			} else {
				is_aromatic = 0;
			}
			//Get H-data
			let h_data = get_h(atom_str, atom_i);
			explH_count = get_h(atom_str, atom_i)['explH_count'];
			//Add data on hydrogens
			if (explH_count > 0) {
				for (let k = 0; k < explH_count; k++) {
					let atom_h = {"atom_id": atom_id + 0.1*(k+1), "atom_symbol": "H", "is_charged": 0, "charge": 0, "atom_start": atom_i, "atom_end": atom_i, "bonds": [{"from_id": atom_id + 0.1*(k+1), "from_pos": atom_start, "to_pos": h_data['position']-1, "type": '-'}] };
					mol_structure_raw.push(atom_h);
				}
			}
			atom_i = h_data['position'];
			//Get charge-data
			charge_data = get_charge(atom_str, atom_i);
			atom.is_charged = charge_data.is_charged;
			atom.charge_qual = charge_data.charge_qual;
			atom.charge_quant = charge_data.charge_quant;
			atom_i = charge_data['position'];
			//Check if it is over
			if (atom_str[atom_i] != ']') {
				alert('Problem while parsing: ' + '\r\n' + curr_section + '\r\n' + 'Please, make corrections or consider alternative input');
				break smiles_traversal;
			}
			//Increase character counter
			i = i + atom_i;
			atom_end = i;
			atom.atom_start = atom_start;
			atom.atom_end = atom_end;
		}
		//Check if current position belongs to the organic atom and characterize it
		if ( bracketed === 0 && organic_set.has(smiles[i]) ) {
			console.log("Organic:  " + i);
			organic = 1;
			//Initialize the variables to store the characteristics of the atom
			atom = {};
			// Basic descriptors
			explH = 0;
			explH_count = 0;
			is_charged = 0;
			is_aromatic = 0;
			charge_qual = "";
			charge_quant = 0;
			is_isotope = 0;
			isotope = "";
			atom_symbol = "";
			atom_bonds = [];
			// Store connections
			let connected_positions = [];
			//positions
			atom_start = i;
			atom_end = i;
			//New atom
			atom_id++;
			atom.atom_id = atom_id;
			atom.bonds = [];
			//Get the substr corresponding to this section
			let atom_str = smiles.substring(i, i + 1);
			//Symbol
			atom_symbol = atom_str;
			if (aromatics_set.has(atom_symbol)) {
				is_aromatic = 1;
			} else {
				is_aromatic = 0;
			}
			//Add the obvious
			atom.is_isotope = is_isotope;
			atom.atom_symbol = atom_symbol;
			atom.is_aromatic = is_aromatic;
			atom.is_charged = is_charged;
			atom.charge_qual = charge_qual;
			atom.charge_quant = charge_quant;
			atom.atom_start = atom_start;
			atom.atom_end = atom_end;
		}
		//Backward search with regards to the current atom
		if (bracketed === 1 | organic === 1) {
			let prebond_found = 0;
			if (atom_start > 0) {
				//Check for explicit bonds behind this atom
				let prebond_data = pre_bond_get(i, atom_start, smiles, atom_bonds, bonds_set);
				if (typeof prebond_data != "undefined") {
					prebond_found = prebond_data['found'];
					atom_start = prebond_data['atom_start'];
				}
				//Find the previous atom if ")"
				let prev_closebranch_data = prev_atom_get_main(i, atom_start, smiles);
				y = prev_closebranch_data;
				//Add the rslts
				if (prev_closebranch_data != "not_the_case") {
					atom_start = prev_closebranch_data['atom_start'];
					//Add the data on bond
					if (prebond_found === 0) {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_closebranch_data['connected_position'], "type": '-'});
					} else {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_closebranch_data['connected_position'], "type": prebond_data['bond']});
					}
				}
				//Find the previous atom if "("
				let prev_openbranch_data = prev_atom_get_frombranch(i, atom_start, smiles);
				//Check and add
				if (prev_openbranch_data != "not_the_case") {
					atom_start = prev_openbranch_data['atom_start'];
					//Add the data on bond
					if (prebond_found === 0) {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_openbranch_data['connected_position'], "type": '-'});
					} else {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_openbranch_data['connected_position'], "type": prebond_data['bond']});
					}
				}
				//Find prev if "DIGIT"
				let prev_digit_data = prev_atom_get_digit (atom_start, smiles);
				//Check and add
				if (prev_digit_data != "not_the_case") {
					atom_start = prev_digit_data['atom_start'];
					//Add the data on bond
					if (prebond_found === 0) {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_digit_data['connected_position'], "type": '-'});
					} else {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_digit_data['connected_position'], "type": prebond_data['bond']});
					}
				}
				//Find prev if exactly atom
				let prev_atom_data = prev_atom_get_atom(atom_start, smiles, connected_positions);
				//Check and add
				if (prev_atom_data != "not_the_case") {
					atom_start = prev_atom_data['atom_start'];
					//Add the data on bond
					if (prebond_found === 0) {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_atom_data['connected_position'], "type": '-'});
					} else {
						atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": prev_atom_data['connected_position'], "type": prebond_data['bond']});
					}
				}
				atom.atom_start = atom_start;
			}
		}
		//Forward search with regards to the bracketed atom.
		//The main search is doing Backward, here is the extension mainly
		//AND rings 
		if (bracketed === 1 | organic === 1) {
			if (atom_end < smiles.length-1) {
				//')' - end of the branch, nothing interesting ahead, just extend the atom
				if (smiles[i+1] == ')') {
					i++
					atom_end++
				}
				//Digit OR '%' -> ring start and the end of the ring should be found
				if ( (/[0-9]/).test(smiles[atom_end+1]) | smiles[atom_end+1] === "%" ) {
					//Get the substr describing belonging to rings
					let ring_str = next_get_rings(i, atom_end, smiles);
					i = ring_str['position'];
					atom_end = ring_str['atom_end'];
					let ring_start = ring_str['ring_start'];
					ring_str = ring_str['ring_str'];
					//Parse ring str into rings
					let these_rings = parse_ring_str (ring_start, ring_str);
					//If some or all rings are not closed yet, they should be closed
					//Thus, check and close
					for (let k = 0; k < these_rings.length; k++) {
						let this_ring = these_rings[k];
						//Check
						ring_action = check_ring(this_ring, closed_rings);
						//Do smth
						if (ring_action === "this ring is closed") {
							//Proceed with the next ring
							continue;
						} else {
							//Find the end of this ring and establish the connection
							let ring_end = end_ring(atom_end, this_ring, smiles);
							//Reportfailure OR Add the data on connection and close the ring officially
							if (ring_end.end_start < 0) {
								alert('Unable to find the end of the ring: ' + '\r\n' + this_ring.str + '\r\n' + 'Please, make corrections or consider alternative input');
								break smiles_traversal;
							} else {
								closed_rings = closed_rings.union(this_ring.pos);
								atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": ring_end.end_start, "type": '-'});
								for (let k = ring_end.end_start; k <= ring_end.end_end; k++) {
									closed_rings.add(k);
								}
							}
						}
					}
					//Temporary solution, correct latter
					i--
					atom.atom_end = atom_end;
				} else {
					//Explicit bond followed by Digit OR '%' -> ring start and the end of the ring should be found
					if ( bonds_set.has(smiles[atom_end+1]) && ((/[0-9]/).test(smiles[atom_end+2]) | smiles[atom_end+2] === "%") ) {
						console.log("Correct atom end:   " + atom_end);
						//Remember the bond
						let ring_bond = smiles[atom_end+1];
						atom_end++
						i++
						//Get the substr describing belonging to rings
						let ring_str = next_get_rings(i, atom_end, smiles);
						i = ring_str['position'];
						atom_end = ring_str['atom_end'];
						let ring_start = ring_str['ring_start'];
						ring_str = ring_str['ring_str'];
						//Parse ring str into rings
						let these_rings = parse_ring_str (ring_start, ring_str);
						//If some or all rings are not closed yet, they should be closed
						//Thus, check and close
						for (let k = 0; k < these_rings.length; k++) {
							let this_ring = these_rings[k];
							//Check
							ring_action = check_ring(this_ring, closed_rings);
							//Do smth
							if (ring_action === "this ring is closed") {
								//Proceed with the next ring
								continue;
							} else {
								//Find the end of this ring and establish the connection
								let ring_end = end_ring(atom_end, this_ring, smiles);
								//Reportfailure OR Add the data on connection and close the ring officially
								if (ring_end.end_start < 0) {
									alert('Unable to find the end of the ring: ' + '\r\n' + this_ring.str + '\r\n' + 'Please, make corrections or consider alternative input');
									break smiles_traversal;
								} else {
									closed_rings = closed_rings.union(this_ring.pos);
									atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": ring_end.end_start, "type": ring_bond});
									for (let k = ring_end.end_start; k <= ring_end.end_end; k++) {
										closed_rings.add(k);
									}
								}
							}
						}
						//Temporary solution, correct latter
						i--
						atom.atom_end = atom_end;
					}
				}
			}
		}
		mol_structure_raw.push(atom);
	}
	// Finalize structure
	let structure_in_progress = mol_structure_raw.filter(element => typeof element != 'undefined');
	// Traverse atoms within mol along with the descriptions obtained and reconstruct the list of bonds as the checkable table
	let bond_table = new Map();
	let priority_bonds = new Set(['=', '#', '$', ':', '\\', '/']);
	// Travesre the structure in progress and introduce some corrections
	for (let i = 0; i < structure_in_progress.length; i++) {
		//Atoms and their positions
		let atom_id = structure_in_progress[i]['atom_id'];
		let atom_symbol = structure_in_progress[i]['atom_symbol'];
		//Traverse the bonds
		atom_bonds = structure_in_progress[i]['bonds'];
		for (let k = 0; k < atom_bonds.length; k++) {
			let destination = atom_bonds[k]['to_pos'];
			let type = atom_bonds[k]['type'];
			//Find the corresponding atom
			let n_dest_atoms = 0;
			for (let j = 0; j < structure_in_progress.length; j++) {
				//Treat H differently
				if (structure_in_progress[j]['atom_symbol']!='H' && structure_in_progress[j]['atom_start'] <= destination && structure_in_progress[j]['atom_end'] >= destination ) {
					let bond_id = [atom_id, structure_in_progress[j]['atom_id']]
								   .sort()
								   .join("-");
					//Check if the bond is known
					if (bond_table.has(bond_id)) {
						//Compare bond types
						let current_record = bond_table.get(bond_id);
						if (!priority_bonds.has(current_record['type']) && priority_bonds.has(type)) {
							bond_table.delete(bond_id);
							bond_table.set(bond_id, {"start_id": bond_id.split('-')[0],
												 "end_id": bond_id.split('-')[1],
												 "start_symbol": structure_in_progress.find((element) => element['atom_id'] === parseFloat(bond_id.split('-')[0]))['atom_symbol'],
												 "end_symbol": structure_in_progress.find((element) => element['atom_id'] === parseFloat(bond_id.split('-')[1]))['atom_symbol'],
												 "type": type});
						} 
					} else {
						//Add bond
						bond_table.set(bond_id, {"start_id": bond_id.split('-')[0],
												 "end_id": bond_id.split('-')[1],
												 "start_symbol": structure_in_progress.find((element) => element['atom_id'] === parseFloat(bond_id.split('-')[0]))['atom_symbol'],
												 "end_symbol": structure_in_progress.find((element) => element['atom_id'] === parseFloat(bond_id.split('-')[1]))['atom_symbol'],
												 "type": type});
					}
				}
			}
		}
	}
	bbbonds = bond_table;
	// Prepare the structure, #red problematic; bond_table should be considered
	const mol_structure = {};
	// Atom storage
	const mol_atoms = [];
	// Bonds that are alrady finalized
	for (let i = 0; i < structure_in_progress.length; i++) {
		let structure_elem = {};
		let elem_poss = Array(structure_in_progress[i]['atom_end'], structure_in_progress[i]['atom_start']).sort();
		let elem_id = structure_in_progress[i]['atom_id'];
		let elem_symbol = structure_in_progress[i]['atom_symbol'];
		let elem_charged = structure_in_progress[i]['is_charged'];
		let elem_isotope = structure_in_progress[i]['is_isotope'];
		let elem_charge_quant = structure_in_progress[i]['charge_quant'];
		if (typeof elem_charge_quant === "undefined") {
			elem_charge_quant = 0;
		}
		let elem_bonds = new Set();
		// Add bonds to this atom
		for (const [key, value] of bond_table.entries()) {
  			let bond_ids = key.split('-');
  			if (bond_ids[0] === elem_id.toString()) {
  				elem_bonds.add(key);
  			}
  			if (bond_ids[1] === elem_id.toString()) {
  				elem_bonds.add(key);
  			}
		}
		// Construct this element of the structure
		structure_elem['atom_id'] = elem_id;
		structure_elem['atom_symbol'] = elem_symbol;
		structure_elem['is_charged'] = elem_charged;
		structure_elem['is_isotope'] = elem_isotope;
		structure_elem['charge_quant'] = elem_charge_quant;
		structure_elem['bonds'] = elem_bonds;
		// Add this element to the structure
		mol_atoms.push(structure_elem);
	}
	mol_structure['atoms'] = mol_atoms;
	mol_structure['bonds'] = bond_table;
	return(mol_structure);
}

// Check number of Carbons and Charge

// Export structure to MOL

</script>

smiles_x = "[CH]CCCCCCC";
result = parse_smiles (smiles_x, aromatics_map, aromatics_set, organic_set, bonds_set, smiles_chars, atom_symbols, bonds_map);