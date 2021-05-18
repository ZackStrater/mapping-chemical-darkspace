from src.structure_to_smiles import convert_to_smiles
import re


'''Converts smiles string to molecular structure'''



class Atom:

    def __init__(self, symbol):
        self.symbol = symbol
        self.bonded_to = []
        self.can_bond = False
        self.heteroatom = False
        if self.symbol not in ["C", "c", "H"]:
            self.heteroatom = True
        self.discovered = False
        self.phantom_bonds = None
        self.phantom_atom = False
        self.fragment = None


class Bond:

    def __init__(self, atom, bond_code):
        self.atom = atom
        self.bond_code = bond_code


class MoleculeStructure:

    def __init__(self):
        self.atom_list = []
        self.fragments_list = []

    def assign_fragment_neighbors(self):
        for fragment in self.fragments_list:
            fragment.assign_fragment_atoms()
        for fragment in self.fragments_list:
            for atom in fragment.atom_list:
                for bond in atom.bonded_to:
                    if bond.atom.fragment != fragment and bond.atom.fragment not in fragment.fragment_bonded_to: # TODO
                        fragment.fragment_bonded_to.append(bond.atom.fragment)

    def H_atoms(self):
        expected_bonds_dict = {
            "N": 3, "O": 2, "P": 3, "Si": 4, "S": 2, "C": 4, "B": 3,
            "b": 3, "c": 4, "n": 3, "o": 2, "p": 3, "s": 2
        }

        bonds_values_dict = {
            1: 1, 2: 2, 3: 3, 4: 1.5, 6: 1, 7: 1, 8: 0,
        }

        hydrogen_atoms = 0

        for atom in self.atom_list:
            if atom.symbol in expected_bonds_dict:  # adding hydrogen atoms
                expected_bonds = expected_bonds_dict[atom.symbol]
                actual_bonds = 0
                for bond in atom.bonded_to:
                    actual_bonds += bonds_values_dict[bond.bond_code]
                if actual_bonds < expected_bonds:
                    hydrogen_atoms += (expected_bonds-actual_bonds)
        return hydrogen_atoms

    def mw(self):
        molecular_weight = 0

        element_weights = {
            "N": 14.0067, "O": 15.9994, "P": 30.9738, "Si": 28.086, "S": 32.06,
            "F": 18.9984, "Cl": 35.453, "Br": 79.904,"I": 126.904, "C": 12.011, "B": 10.81,
            "b": 10.81, "c": 12.011, "n": 14.0067, "o": 15.9994, "p": 30.9738, "s": 32.06, "R": 0
        }

        for atom in self.atom_list:
            molecular_weight += element_weights[atom.symbol]

        molecular_weight += 1.00784*self.H_atoms()

        return molecular_weight

    def formula(self):
        symbol_counts = {
            "C": 0, "H": 0, "B": 0, "Br": 0, "Cl": 0, "F": 0, "I": 0, "N": 0, "O": 0, "P": 0,  "S": 0, "Si": 0,
            "b": 0, "c": 0, "n": 0, "o": 0, "p": 0, "s": 0, "R": 0
        }

        for atom in self.atom_list:
            symbol_counts[atom.symbol] += 1

        symbol_counts["H"] = self.H_atoms()

        lowercase_symbols = []
        for key in symbol_counts:
            if key.islower():
                symbol_counts[key.upper()] += symbol_counts[key]
                lowercase_symbols.append(key)

        for symbol in lowercase_symbols:
            del symbol_counts[symbol]

        formula_string = ""

        for key, value in symbol_counts.items():
            if value > 0:
                formula_string += key + str(value)

        return formula_string.translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))



class Fragment(MoleculeStructure):
    def __init__(self, name, atom_list):
        MoleculeStructure.__init__(self)
        self.name = name
        self.atom_list = atom_list
        self.fragment_bonded_to = []

    def generalize_heterocycle_name(self):
        heteroatoms = 0
        for atom in self.atom_list:
            if atom.heteroatom:
                heteroatoms += 1
        self.name = f'{heteroatoms}-' + self.name

    def assign_fragment_atoms(self):
        for atom in self.atom_list:
            atom.fragment = self


bond_encoder = {
    "=": 2,  # double bond
    "#": 3,  # triple bond
    "%": 4, # aromatic bonds also denoted by lowercase element symbols
    "$": 5,  # ionic bond
    "/": 6,  # vinyl up
    "\\": 7,  # vinyl down
    ".": 8,  # non-bond
    "&": 9  # any bond
}


def encode_bond(bonding_info):
    code = None
    for key in bond_encoder:
        if key in bonding_info:
            code = (bond_encoder[key])
    if code is None:
        return 1
    else:
        return code


def convert_to_structure(molecule, smiles_string):

    pre_correction = re.sub(r"|\(H\)|\(\[H\]\)", "", smiles_string)
    corrected_smiles_string = re.sub(r"[\[\]H@+-]", "", pre_correction)
    # if fragment has phantom atoms (signified by portions encapsulated by {}), add an * to each of those atoms
    if re.search(r"[{}]", corrected_smiles_string):
        def add_phantom_atoms(matchobj):
            first_edit = re.sub(r"(N|O|P|Si|S|F|Cl|Br|I|C|B|R|Q|b|c|n|o|p|s)", r"\1*", matchobj.group())
            return re.sub(r"[{}]", r"", first_edit)

        corrected_smiles_string = re.sub(r"({)(\S*?)(})", add_phantom_atoms, corrected_smiles_string)
    for match in re.findall(r"N|O|P|Si|S|F|Cl|Br|I|C|B|E|b|c|n|o|p|s|R|Q|W|X|Y|Z", corrected_smiles_string):
        molecule.atom_list.append(Atom(match))
        # create Atom object for each elemental symbol found in the smiles_string, keeps order of Atom in the string
        # R = any element
        # Q = any heteratom
    bond_map = re.findall(r"(?:N|O|P|Si|S|F|Cl|Br|I|C|B|E|b|c|n|o|p|s|R|Q|W|X|Y|Z)([^A-Za-z]*)", corrected_smiles_string)
    # ordered list containing all bonding symbol denoting bonding information following each element symbol

    left_parens_list = []
    # ordered list of atoms that are followed by a "("
    ring_closure_dict = {}
    # keeps track of ring closures

    for i in range(0, len(molecule.atom_list)):
        # for loop through all Atoms, but allows access to Atom (i + 1)
        # notice this range includes the last item in the list cause ring closure can be included in last atom
        all_bond_info = (bond_map[i])
        # string with all bonding info for current atom

        # atoms followed by "*" in the string are considered phantom atoms - don't get marked as discovered and
        # can be found in areas of the molecule that have already been assigned to a fragment
        if "*" in all_bond_info:
            molecule.atom_list[i].phantom_atom = True

        ring_closure_info = re.findall(r"(\D*[%]\d{2}|\D*\d)", all_bond_info)
        # pulls out all ring #'s and the bond info that proceeds them

        for ele in ring_closure_info:
            ring_closure_split = re.split(r"(\D*(?=\d))", ele, 1)
            # separates bonding info [1] from ring # [2]
            if ring_closure_split[2] in ring_closure_dict:
                # if this ring # has already been encountered -> make bond between those two Atoms
                closure_partner = ring_closure_dict[ring_closure_split[2]]
                combined_bonding_info = closure_partner[1] + ring_closure_split[1]
                # often bonding info only included for one of the ring closure partners, so combine and use for both Bonds
                molecule.atom_list[i].bonded_to.append(Bond(closure_partner[0], encode_bond(combined_bonding_info)))
                closure_partner[0].bonded_to.append(Bond(molecule.atom_list[i], encode_bond(combined_bonding_info)))
                # adding bonds to both closure partners
            else:
                ring_closure_dict[ring_closure_split[2]] = (molecule.atom_list[i], ring_closure_split[1])
                # enter atom in ring_closure_dict -> ring #: (Atom, bonding info)

        bond_info = re.split(r"(\D*$)", all_bond_info, 1)[1]
        # captures the rest of bonding info after ring closure #'s
        if ")" not in bond_info:
            # make Bond between current Atom and Atom (i + 1)
            if i != (len(molecule.atom_list) - 1):
                # prevents index error on last Atom
                molecule.atom_list[i].bonded_to.append(Bond(molecule.atom_list[i + 1], encode_bond(bond_info)))
                molecule.atom_list[i + 1].bonded_to.append(Bond(molecule.atom_list[i], encode_bond(bond_info)))
                if "(" in bond_info:
                    left_parens_list.append(molecule.atom_list[i])

        else:
            # if ")" in bond_info, count the # of ")"'s and go back that many "("'s in the string
            # and make Bond from atom before last "(" and the next atom (i + 1)
            if i != (len(molecule.atom_list) - 1):
                # important not to include last item here, ")" at the end of a string don't actually give structural information
                parentheses_count = bond_info.count(")")
                # counts # of ")" in bond_info and bonds current Atom to
                bond_atom1 = left_parens_list[-parentheses_count]
                # atom that opened the original bracket
                bond_atom2 = molecule.atom_list[i + 1]
                # atom that comes after all the ")"'s
                bond_atom2.bonded_to.append(Bond(bond_atom1, encode_bond(bond_info)))
                bond_atom1.bonded_to.append(Bond(bond_atom2, encode_bond(bond_info)))
                if "(" in bond_info:  # bond info has "C)(C"
                    for c in range(parentheses_count - 1):
                        left_parens_list.pop()
                        # this keeps track of last add "C(" instance
                        # i.e. A(B(C))(DE)F , for bond info between C and D -> C))(D
                        # we get two ")" so we should delete two atoms off left_parens_list
                        # but we also have "(" so we need to add an atom back on the list
                        # the correct atom is not the previous one, but actually the second one we deleted from left_parens_list
                        # so instead we just delete one less from the list and now it will correctly bond A to F

                else:
                    for c in range(parentheses_count):
                        left_parens_list.pop()

    phantom_bonds_dict = {"W": 1, "X": 2, "Y": 3, "Z": 4}
    # used to decipher number of phantom bonds on atom

    for atom in molecule.atom_list:
        if atom.symbol.islower():
            for bond in atom.bonded_to:
                if bond.atom.symbol.islower():
                    bond.bond_code = 4
    # smiles strings sometimes use lowercase to denote an aromatic ring
    # without explicit bond notation, encode_bond will return a 1 (single bond) for bond.code
    # this loop records the aromatic bonding information for strings with such notation,
    # only counting bonds between two lower case atom symbols

        for bond in atom.bonded_to:
            if bond.atom.symbol == "E":
                atom.can_bond = True
                # here the E atoms represent sites that the fragment can bond
                # first atom is bonded to element E, then can_bond for atom changed to True, then delete [R] atoms from
                # molecule.atom_list and any bonds to E atoms

            if bond.atom.symbol in phantom_bonds_dict:
                atom.phantom_bonds = phantom_bonds_dict[bond.atom.symbol]
            # phantom bonds allows checking number of bonds atom should have without traversing to those atoms
            # for distinguishing between carboxylic acid and ester, or between primary, secondary, tertiary amines
        atom.bonded_to = [bond for bond in atom.bonded_to if bond.atom.symbol != "E" and bond.atom.symbol not in phantom_bonds_dict]
    molecule.atom_list = [atom for atom in molecule.atom_list if atom.symbol != "E" and atom.symbol not in phantom_bonds_dict]
    return molecule



