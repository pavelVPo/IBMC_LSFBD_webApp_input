// 16-11-2025, the problem with the percieving of aromaticity has been solved,
// testing is requaired:
// Cn1cc(-c2cc3c(N4CCN(c5ncc([C@@](C)(N)c6ccc(F)cc6)cn5)CC4)ncnn3c2)cn1
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


// 1	Remove stereo and associated hydrogens, since PASS does not consider it
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
	smiles = smiles.replaceAll("[B]", "B");
	smiles = smiles.replaceAll("[C]", "C");
	smiles = smiles.replaceAll("[N]", "N");
	smiles = smiles.replaceAll("[O]", "O");
	smiles = smiles.replaceAll("[P]", "P");
	smiles = smiles.replaceAll("[S]", "S");
	smiles = smiles.replaceAll("[F]", "F");
	smiles = smiles.replaceAll("[Cl]", "Cl");
	smiles = smiles.replaceAll("[Br]", "Br");
	smiles = smiles.replaceAll("[I]", "I");
	smiles = smiles.replaceAll("[CH]", "C");
	smiles = smiles.replaceAll("[=CH]", "=C");
	smiles = smiles.replaceAll("[#CH]", "#C");
	smiles = smiles.replaceAll("[$CH]", "$C");
	smiles = smiles.replaceAll("[:CH]", ":C");
	smiles = smiles.replaceAll("[\\CH]", "\\C");
	smiles = smiles.replaceAll("[/CH]", "/C");
	smiles = smiles.replaceAll("[cH]", "c");
	smiles = smiles.replaceAll("[:cH]", "c");
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
			//Count '(' and ')'
			if (smiles[k] === ')') {
				open_branch++;
			}
			if (smiles[k] === '(') {
				open_branch--;
			}
			if (open_branch === 0 && smiles[k-1] != ')') {
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
				if (open_branch === 0 && smiles[k-1] != ')') {
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

// 13 Check if this ring is over
function check_ring (this_ring, closed_rings) {
	if ( closed_rings.intersection(this_ring.pos).size === 0 ) {
		return "close this ring";
	} else {
		return "this ring is closed";
	}
}

// 14 Get the ending positon of this ring
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

// 15 Get the starting position of the atom starting the ring given the position of the first numeric or % describing the ring
// for aromatics only
function extend_start (ring_pseudo_start, smiles) {
	let ring_start = ring_pseudo_start;
	// Just decrease the counter until the aromatic atom is found
	for (let i = ring_pseudo_start; i > 0; i--) {
		if ( aromatics_set.has(smiles[i-1] + smiles[i-2]) ) {
			ring_start = i-2;
			return ring_start;
		} else {
			if ( aromatics_set.has(smiles[i-1]) ) {
				ring_start = i-1;
				return ring_start;
			}
		}
	}
	alert("Smth is wrong");
	return ring_pseudo_start;
}

// 16 Get all the rings from SMILES
function get_all_aromatic_rings (smiles) {
	// Initialize some variables
	const closed_rings = new Set();
	let rings = [];
	let position = 0;
	// Traverse the smiles detecting rings' borders
	const ring_marks = parse_ring_str(0, smiles);
	// Collect the corresponding borders to records
	for (let i = 0; i < ring_marks.length; i++) {
		const ring_id = ring_marks[i]['str'];
		let ring_start = [...ring_marks[i]['pos']].toSorted()[0];
		ring_start = extend_start(ring_start, smiles);
		const ring_start_fin = [...ring_marks[i]['pos']].toSorted()[[...ring_marks[i]['pos']].length-1];
		// Find the end of this ring if it is not ended still
		if (!closed_rings.has(ring_id)) {
			for (let k = 0; k < ring_marks.length; k++) {
				const ring_id_end = ring_marks[k]['str'];
				const ring_end = [...ring_marks[k]['pos']].toSorted()[[...ring_marks[k]['pos']].length-1];
				// Pair it
				if (ring_id === ring_id_end && i != k && ring_start < ring_end) {
					let this_ring = {"id": i.toString() + "_" + ring_id, "start": ring_start, "end": ring_end}
					rings.push(this_ring);
				}
			}
		}
	}
	return rings;
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
	let point_dlm;
	point_dlm = input_smiles.indexOf('.');
	if (point_dlm > -1) {
		let smiles_array = input_smiles.split('.');
		input_smiles = smiles_array.reduce( function (a, b) { return a.length > b.length ? a : b; } ).trim();
	}
	// Check the length
	if (input_smiles.length > 500) {
		return("smth_wrong");
		alert("Input string is too long");
	}
	if (input_smiles.length < 3) {
		alert("Input string is too short");
		return("smth_wrong");
	}
	// Return the SMILES string
	return(input_smiles);
}

////////////////////////////////////////////////////////////////////////////////
// Functions to process SMILES												///
//////////////////////////////////////////////////////////////////////////////

// 1. Remove flags on stereochemistry, and check length again
function destereo_smiles (smiles_input) {
	const smiles = clear_stereo(smiles_input);
	if (smiles.length < 3) {
		alert("Input SMILES string is too short");
		return("smth_wrong");
	}
	return(smiles);
}

// 2. Parse SMILES,
// #red this function is still immature and should be checked latter: some unused vars could be present and some exceptions could be unsecured, some variables have questionable names, some loging is no longer needed
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
			atom.is_aromatic = is_aromatic;
			//Get H-data
			let h_data = get_h(atom_str, atom_i);
			explH_count = get_h(atom_str, atom_i)['explH_count'];
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
			//Add data on hydrogens
			atom_id++
			if (explH_count > 0) {
				for (let k = 0; k < explH_count; k++) {
					atom_id++
					let atom_h = {"atom_id": atom_id, "atom_symbol": "H", "is_charged": 0, "is_aromatic": 0, "charge_quant": 0, "charge_qual": "", "atom_start": i - atom_i + h_data['position'], "atom_end": i - atom_i + h_data['position'], "bonds": [{"from_id": atom_id, "from_pos": atom_start, "to_pos": i - atom_i + h_data['position'] - 1, "type": '-'}] };
					mol_structure_raw.push(atom_h);
				}
			}
		}
		//Check if current position belongs to the organic atom and characterize it
		if ( bracketed === 0 && (organic_set.has(smiles[i]) | organic_set.has(smiles[i] + smiles[i+1])) ) {
			organic = 1;
			let atom_str;
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
			if ( organic_set.has(smiles[i] + smiles[i+1]) ) {
				atom_end = i + 1;
				//Get the substr corresponding to this section
				atom_str = smiles.substring(i, i + 2);
			} else {
				atom_end = i;
				//Get the substr corresponding to this section
				atom_str = smiles.substring(i, i + 1);
			}
			//New atom
			atom_id++;
			atom.atom_id = atom_id;
			atom.bonds = [];
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
					let these_rings = parse_ring_str(ring_start, ring_str);
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
								//Check if the number marking the end of this ring is preceeded by the explicit bond
								if (bonds_set.has(smiles[ring_end.end_start-1])) {
									atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": ring_end.end_start, "type": smiles[ring_end.end_start-1]});
								} else {
									atom.bonds.push({"from_id": atom_id, "from_pos": atom_start, "to_pos": ring_end.end_start, "type": '-'});
								}
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
	sip = structure_in_progress;
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
												 "start_symbol": structure_in_progress.find((element) => element['atom_id'] === parseInt(bond_id.split('-')[0]))['atom_symbol'],
												 "end_symbol": structure_in_progress.find((element) => element['atom_id'] === parseInt(bond_id.split('-')[1]))['atom_symbol'],
												 "type": type});
						} 
					} else {
						//Add bond
						bond_table.set(bond_id, {"start_id": bond_id.split('-')[0],
												 "end_id": bond_id.split('-')[1],
												 "start_symbol": structure_in_progress.find((element) => element['atom_id'] === parseInt(bond_id.split('-')[0]))['atom_symbol'],
												 "end_symbol": structure_in_progress.find((element) => element['atom_id'] === parseInt(bond_id.split('-')[1]))['atom_symbol'],
												 "type": type});
					}
				}
			}
		}
	}
	// Prepare the structure
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
		let elem_charge_qual = structure_in_progress[i]['charge_qual'];
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
		structure_elem['charge_qual'] = elem_charge_qual;
		structure_elem['bonds'] = elem_bonds;
		// Add this element to the structure
		mol_atoms.push(structure_elem);
	}
	mol_structure['atoms'] = mol_atoms;
	mol_structure['bonds'] = bond_table;
	return(mol_structure);
}

// 3. Check the number of Carbons
function check_3c (mol_structure) {
	let condition_3c = 3;
	for (let i = 0; i < mol_structure['atoms'].length; i++) {
		if (mol_structure['atoms'][i]['atom_symbol'] === "C" | mol_structure['atoms'][i]['atom_symbol'] === "c") {
			condition_3c--;
		}
		if (condition_3c <= 0) {
			return true;
		}
	}
	return false;
}

// 4. Reconstruct rings, it could be done and probably should be done, during the initial SMILES parsing, but since
// performance is not the main concern at the moment
// function to parse SMILES is quite big and will be reworked latter
// let it be the distinct function
/*function reconstruct_aromatics_test (mol_structure, all_rings, smiles) {
	// Count aromatic rings
	ring_number = 0;
	// Atom in rings should have at least two bonds and be aromatic
	for (let value of mol_structure['atoms']) {
		if ( aromatics_set.has(value['atom_symbol']) && value['bonds'].size > 1 ) {
			// The main code goes here
			// Find the ring for this atom:
			// - this atom could be the start OR the end of this ring

		} else {
			if ( aromatics_set.has(value['atom_symbol']) && value['bonds'].size < 2 ) { 
				return("smth_wrong");
				alert("Aromatic atoms not in rings are not allowed");	
			}
		}
	}
}*/

// 5. Reconstruct aromatic bonds and mark those, which are in rings and not between them
// Bond is in ring if
// its atoms belong to single ring
// its atoms belong to two rings
function reconstruct_aromatics (mol_structure) {
	const aromatic_adjlist = new Map();
	for (let i = 0; i < mol_structure['atoms'].length; i++) {
		const this_atom_pairs = new Set();
		if ( aromatics_set.has(mol_structure['atoms'][i]['atom_symbol']) ) {
			for (let value of mol_structure['atoms'][i]['bonds']) {
				let this_pair = parseInt(value.split("-").filter((element) => parseInt(element) != mol_structure['atoms'][i]['atom_id']));
				// Filter to get corrsponding array element and atom symbol from it
				let this_pair_symbol = mol_structure['atoms'].filter( (element) => element['atom_id'] === this_pair )[0]['atom_symbol'];
				if (aromatics_set.has(this_pair_symbol)) {
					this_atom_pairs.add(this_pair);
				}
			}
			if (mol_structure['atoms'][i]['atom_id'] === 2) {
			}
			aromatic_adjlist.set(mol_structure['atoms'][i]['atom_id'], this_atom_pairs);
		}
	}
	x = aromatic_adjlist;
	for (let [key, value] of mol_structure['bonds']) {
		let aromatic_status = value;
		// Make all bonds non-aromatic 
		aromatic_status['aromatic'] = 0;
		mol_structure['bonds'].set(key, aromatic_status);
		if (aromatics_set.has(value['start_symbol']) && aromatics_set.has(value['end_symbol'])) {
			// Check if this bond is in aromatic ring, i.e. without this bond its atoms are still connected via aromatics
			// Breadth-first search ignoring the immediate connection to the target node (atom)
			let queue = [];
			const observed_atoms = new Set();
			// Initialize search ignoring the target atom
			queue.push(parseInt(value['start_id']));
			observed_atoms.add(parseInt(value['start_id']));
			// Actual Breadth-first search
			let preventive_counter = 0;
			do {
				preventive_counter++;
				console.log("search:");
				console.log(value['start_id']);
				console.log(value['end_id']);
				console.log(queue);
				let current_start = queue.shift();
				console.log(current_start);
				// Check if this the target atom
				if (current_start === parseInt(value['end_id'])) {
					aromatic_status['aromatic'] = 1;
					mol_structure['bonds'].set(key, aromatic_status);
					break;
				}
				// Proceed
				// Copying from the Set is a big deal
				let adj_toSearch = new Set([...aromatic_adjlist.get(current_start)]);
				// Delete the target atom on the first step
				if (preventive_counter === 1) {
					adj_toSearch.delete(parseInt(value['end_id']));
				}
				for (let setval of adj_toSearch) {
					if (!observed_atoms.has(setval)) {
						queue.push(setval);
						observed_atoms.add(setval);
					}
				}
			} while (queue.length > 0);
		} 
	}
	return mol_structure;
}

// Export structure to MOL just for PASS, coordinates will not be computed
function structure_to_mol (mol_structure) {
	// What is MOL? It is the format to store chemical structures; which PASS could read
	// SEE: https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf
	// MOL has more or lesse arbitrary header block spanning accross 3 lines  
	const header_block = "\r\nJust a basic MOL for PASS only from JS without chirality and all\r\n";
	// MOL has counts line with number of atoms, etc. To fill this line it is needed to be known number of lines for additional properties.
	// The format is aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
	// aaa -> number of atoms
	// bbb -> number of bonds
	// lll -> number of atom lists
	// fff -> obsolete
	// ccc -> chiral flag; 0=not chiral, 1=chiral
	// sss -> obsolete
	// xxx -> obsolete
	// rrr -> obsolete
	// ppp -> obsolete
	// iii -> obsolete
	// mmm -> 999
	const number_of_atoms = mol_structure['atoms'].length;
	const number_of_bonds = mol_structure['bonds'].size;
	const zero_trail = "   ";
	const zero_val = "  0";
	const counts_line = [zero_trail,number_of_atoms].join("").slice(-3) + [zero_trail,number_of_bonds].join("").slice(-3) + zero_val + zero_val + zero_val + zero_val + zero_val + zero_val + zero_val + zero_val + "999" + " V2000";
	// MOL has the atom block, Atom Block is made up of atom lines, one line per atom
	// The format is xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
	// xyz -> atom coordinates
	// aaa -> atom symbol
	// dd -> mass difference (0)
	// ccc -> charge, very interesting: (0 if uncharged; 1, 2, 3 if positive charges where 1= +3, 2 = +2, 3 = +1; 4 if doublet radical; 5, 6, 7 if negative charges where 5 = -1, 6 = -2, 7 = -3)
	// sss -> atom stereo parity, 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked stereo center
	// hhh -> hydrogen count + 1; H0 - no hydrogens unless explicitly drawn
	// bbb -> ignore stereo care
	// vvv -> valence; 0 - no marking
	// HHH -> H0 designator; 0 - not specified
	// rrr, not used
	// iii, not used
	// mmm -> for reactions
	// nnn -> for reactions
	// eee -> for rezctions
	const atom_block = [];
	// MOL has the bond block, The Bond Block is made up of bond lines, one line per bond, with the following format:
	// 111222tttsssxxxrrrccc
	// 111 first atom number 1 - number of atoms
	// 222 second atom number 1 - number of atoms
	// ttt bond type 1 = Single, 2 = Double, 3 = Triple, 4 = Aromatic, 5 = Single or Double, 6 = Single or Aromatic, 7 = Double or Aromatic, 8 = Any
	// sss bond stereo Single bonds: 0 = not stereo, 1 = Up, 4 = Either, 6 = Down, Double bonds: 0 = Use x-, y-, z-coords from atom block to determine cis or trans, 3 = Cis or trans (either) double bond
	// xxx not used
	// rrr bond topology 0 = Either, 1 = Ring, 2 = Chain
	// ccc -> reactions only
	bond_block = [];
	// MOL has the properties block;
	// The Properties Block is made up of mmm lines of additional properties, Note: mmm is no longer supported and is set to 999 as the default.
	// The prefix: M END terminates the properties block.
	// Charge, it superseeds all charge and radical value
	// Radical, his property supersedes all charge and radical values in the atom
	// Isotope,
	// Ring Bond Count, Number of ring bonds allowed
	// Substitution Count
	// Unsaturated Count
	// Link atom
	// ...
	// M END -> end of block
	let n_of_charged = 0;
	let charge_line = "";
	const prop_block = [];
	// So, to construct MOL it is needed to go through the structure describe and count atom and bonds
	// Filling the atom_block, bond_block and prop_block
	// Since bonds are between the pairs of atoms, it is necessary to monitor them
	described_bonds = new Set();
	for (let i = 0; i < mol_structure['atoms'].length; i++) {
		// Construct MOL
		// Construct this atom's line
		let atom_line;
		// Populate this atom's line in atom block
		// The format is xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
		let atom_coords = "    0.0  00   0.0  0    0.0  0 ";
		let atom_symbol = mol_structure['atoms'][i]['atom_symbol'];
		// convert the atom symbol, IF aromatic
		if (aromatics_set.has(atom_symbol)) {
			if ( atom_symbol === 'c') {atom_symbol = 'C'} else {
				if ( atom_symbol === 'b') {atom_symbol = 'B'} else {
					if ( atom_symbol === 'n') {atom_symbol = 'N'} else {
						if ( atom_symbol === 'o') {atom_symbol = 'O'} else {
							if ( atom_symbol === 'p') {atom_symbol = 'P'} else {
								if ( atom_symbol === 's') {atom_symbol = 'S'} else {
									if ( atom_symbol === 'as') {atom_symbol = 'As'} else {
										if ( atom_symbol === 'se') {atom_symbol = 'Te'} else {
											if ( atom_symbol === 'se') {atom_symbol = 'Se'} else {
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		let atom_massdif = "   0";
		let atom_charge = "  0";
		let atom_stereo = "  0";
		let atom_h = "  0";
		let atom_b = "  0";
		let atom_v = "  0";
		let atom_H = "  0";
		let atom_r = "  0";
		let atom_i = "  0";
		let atom_m = "  0";
		let atom_n = "  0";
		let atom_e = "  0";
		// Gather to the atom line
		atom_line = atom_coords + atom_symbol + atom_massdif + atom_charge + atom_stereo + atom_h + atom_b + atom_v + atom_H + atom_r + atom_i + atom_m + atom_n + atom_e;
		// Add this line to this atom's block
		atom_block.push(atom_line);
		
		// Construct bond lines for this atom considering the bonds, which exist
		mol_structure['atoms'][i]['bonds'].forEach((bond_i_key, bond_i) => {
			if (!described_bonds.has(bond_i_key)) {
				let bond_line;
				// Populate this atom's line in bond block
				// The format is 111222tttsssxxxrrrccc
				// Consider hydrogens
				let atom_one_number = mol_structure['atoms'].findIndex((element) => element['atom_id'] === parseInt(bond_i_key.split("-")[0])) + 1;
				let atom_two_number = mol_structure['atoms'].findIndex((element) => element['atom_id'] === parseInt(bond_i_key.split("-")[1])) + 1;
				atom_one_number = [zero_trail, atom_one_number].join("").slice(-3);
				atom_two_number = [zero_trail, atom_two_number].join("").slice(-3);
				let bond_type = mol_structure['bonds'].get(bond_i)['type'];
				// The codes for bond types are: 1 = Single, 2 = Double, 3 = Triple, 4 = Aromatic, 5 = Single or Double, 6 = Single or Aromatic, 7 = Double or Aromatic, 8 = Any
				// Convert the bond type to the format appropriate for MOL
				// Consider aromatics, #red this is not really handy, the type must be specified explicitly
				if (bond_type === "-" && mol_structure['bonds'].get(bond_i)['aromatic'] === 1) {
					bond_type = "  4";
				} else {
					if (bond_type === "-" | bond_type === "\\" | bond_type === "/") {
						bond_type = "  1";
					} else {
						if (bond_type === "=") {
							bond_type = "  2";
						} else {
							if (bond_type === "#") {
								bond_type = "  3";
							} else {
								if (bond_type === ":") {
									bond_type = "  4";
								}
							}
						}
					}
				}
				let bond_stereo = "  0";
				let bond_x = "   ";
				let bond_r = "  0";
				let bond_c = "  0";
				// Gather to the bond line
				bond_line = atom_one_number + atom_two_number + bond_type + bond_stereo + bond_x + bond_r + bond_c;
				bond_block.push(bond_line);
				described_bonds.add(bond_i_key);
			}
		});

		// Construct property lines for this atom, charge is an important property in this case
		// Format in the most cases is M  XXX
		// So, if atom is charged -> add the property line, LIKE:
		// M__CHG  1____8__-1: prefix, number of entries on the line, charge [-15, +15]
		// #red add the check-up for charge belonging to [15, 15]
		let atom_charge_q = mol_structure['atoms'][i]['charge_qual'] + mol_structure['atoms'][i]['charge_quant'].toString();
		let atom_number_here = ["    ", (i+1).toString()].join("").slice(-4);
		if (atom_charge_q != "0" ) {
			n_of_charged++
			// Stringify the charge value
			atom_charge_q = ["    ", atom_charge_q].join("").slice(-4);
			// Construct property line
			charge_line = charge_line + atom_number_here + atom_charge_q;
		} 

	}
	charge_line = "M  CHG" + ["    ", n_of_charged.toString()].join("").slice(-3) + charge_line;
	if (n_of_charged > 0) {
		prop_block.push(charge_line);
		prop_block.push("\r\n");
	}
	bond_block = [... new Set(bond_block)];
	// Construct MOL string
	let mol_string = header_block + "\r\n" + counts_line + "\r\n" + atom_block.join("\r\n") + "\r\n" + bond_block.join("\r\n") + "\r\n"+ prop_block.join("") + "M  END\r\n";
	return mol_string;
}