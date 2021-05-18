import re

'''Converts molecular structure to smiles string'''


bond_decoder = {
    1: "",  # single bond
    2: "=",  # double bond
    3: "#",  # triple bond
    4: "",  # aromatic bond
    5: "$",  # ionic bond
    6: "/",  # vinyl up
    7: "\\",  # vinyl down
    8: ".",  # non-bond
    9: "&" # any bond
    }


def decode_bond(code):
    return bond_decoder[code]


class AtomSmiles:
    # data object for atoms.  keeps track of bonding information, ring closure bond info, and closure partners

    def __init__(self):
        self.bond_code = None
        self.ring_closures = []
        self.closure_partners = []


def convert_to_smiles(molecule, start):
    smiles_string = ""
    smiles_construction_list = []
    # ordered list of atoms and parentheses
    atom_info = {}
    # dictionary linking atoms to their AtomSmiles info
    visited = {}
    # keeping track of which atoms have already been visited
    ring_counter = 1
    # keeps track and numbers rings
    for a in molecule.atom_list:
        visited[a] = False
    # all atoms start unvisited

    def dfs(current_atom, previous_atom):
        nonlocal ring_counter
        visited[current_atom] = True
        smiles_construction_list.append(current_atom)

        previous_bond = None
        for bond in current_atom.bonded_to:
            if bond.atom == previous_atom:
                previous_bond = bond
        # this algorithm is backward looking, i.e. bond info

        atom_info[current_atom] = AtomSmiles()
        if previous_bond:
            atom_info[current_atom].bond_code = previous_bond.bond_code

        # ring detection
        for bond in current_atom.bonded_to:
            if bond != previous_bond and visited[bond.atom] and bond.atom not in atom_info[current_atom].closure_partners:
                # ring detected if the bond was not the bond we just came from
                # atom has been visited already
                # and atom has not previously marked as a closure partner with current atom (prevents double counting rings)
                atom_info[current_atom].ring_closures.append(decode_bond(bond.bond_code))
                if ring_counter > 9:
                    atom_info[current_atom].ring_closures.append("%")
                    atom_info[bond.atom].ring_closures.append("%")
                atom_info[current_atom].ring_closures.append(ring_counter)
                atom_info[bond.atom].ring_closures.append(ring_counter)
                # formatting (bond info)(% if ring # > 9)(ring #)
                # bonding info only applied to current atom, not both closure partners (accepted smiles notation)
                atom_info[current_atom].closure_partners.append(bond.atom)
                atom_info[bond.atom].closure_partners.append(current_atom)

                ring_counter += 1

            untraveled_paths = []
            for b in current_atom.bonded_to:
                if not visited[b.atom]:
                    untraveled_paths.append(bond.atom)
            if not visited[bond.atom]:
                if len(untraveled_paths) > 1:
                    smiles_construction_list.append("(")
                dfs(bond.atom, current_atom)
                if len(untraveled_paths) > 1:
                    smiles_construction_list.append(")")
                # parentheses added to denote fragments attached to current atom
                # parentheses not needed if only one path is left

    dfs(start, None)
    # begin dfs process

    # atom_ring_list = []
    # for ele in smiles_construction_list:
    #     if isinstance(ele, Atom):
    #         if len(atom_info[ele].ring_closures) > 0:
    #             atom_ring_list.append(smiles_construction_list.index(ele))

    # building the smiles string
    for ele in smiles_construction_list:
        if isinstance(ele, str):
            smiles_string += ele
        else:
            if atom_info[ele].bond_code:
                smiles_string += decode_bond(atom_info[ele].bond_code)
            smiles_string += ele.symbol
            for c in atom_info[ele].ring_closures:
                smiles_string += str(c)

    return smiles_string






