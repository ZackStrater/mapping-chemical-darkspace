from smiles_to_structure import convert_to_structure, MoleculeStructure
# from fragments_building_library import heterocycles
from collections import Counter


class AtomData:
    def __init__(self, symbol):
        self.symbol = symbol
        self.bonds = []
        self.daughter_branches = []


class Branch:
    def __init__(self, branch_num):
        self.branch_num = branch_num
        self.sequence = []


class FragmentMap:
    def __init__(self):
        self.branches = {}
        self.base_atom = None

    def show_map(self):
        print("base atom:")
        print(self.base_atom.symbol)
        print("base_atom branches:")
        for branch in self.base_atom.daughter_branches:
            print("branch number:")

            print(branch[0])


def abbr_bond(bond):
    return bond.bond_code, bond.atom.symbol
# tuple representation of a bond


def fragmentize(smiles_string, fragment_string):  #specify fragment list

    total_molecule = convert_to_structure(MoleculeStructure(), smiles_string)
    current_fragment = convert_to_structure(MoleculeStructure(), fragment_string)

    def find_atom(fragment):
        for ele in ["Si", "P", "I", "Br", "Cl", "F", "B", "S", "O", "N", "C"]:
            for atom in fragment.atom_list:
                if atom.symbol == ele:
                    return atom
    base_atom = find_atom(current_fragment)

    def map_fragment(fragment, base_atom):
        fragment_map = FragmentMap()
        # keeps track of branch connectivity
        visited = {}
        for atom in fragment.atom_list:
            visited[atom] = False
        # keeps track of which atoms have been visited
        branch_num = 1  # TODO this might be unnecessary with daughter_branches

        def dfs(current_atom, previous_atom, branch):
            visited[current_atom] = True

            current_atom_data = AtomData(current_atom.symbol)
            # data object for current atom

            if branch:
                branch.sequence.append(current_atom_data)
            else:
                fragment_map.base_atom = current_atom_data
            # append atom data to current branch
            # if not branch, current atom is base_atom, add as beginning part of map

            previous_bond = None
            for bond in current_atom.bonded_to:
                if bond.atom == previous_atom:
                    previous_bond = bond

            if len(current_atom.bonded_to) > 2:
                # data objects for atoms added sequentially to branch.sequence
                for bond in current_atom.bonded_to:
                    if bond != previous_bond and not visited[bond.atom]:
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)
                        # current_atom will be the previous atom for new branch
            # if atom is a branch point, create new branch for each path

            elif len(current_atom.bonded_to) == 2:
                for bond in current_atom.bonded_to:
                    if bond != previous_bond:
                        current_atom_data.bonds.append(abbr_bond(bond))

                for bond in current_atom.bonded_to:
                    if not visited[bond.atom]:
                        dfs(bond.atom, current_atom, branch)
            # elif atom is contiguous point, continue to next atom

            else:
                if not previous_bond:
                    for bond in current_atom.bonded_to:
                        current_atom_data.bonds.append(abbr_bond(bond))
                    for bond in current_atom.bonded_to:
                        if not visited[bond.atom]:
                            dfs(bond.atom, current_atom, branch)
                # else atom is either an endpoint, in which case branch ends
                # if atom does not have previous bond, it means that mapping started from an atom with only 1 bond

        def new_branch(current_atom, previous_atom, previous_atom_data, bond_code):
            # current atom used to start dfs chain
            # previous atom helps tell dfs where the previous atom was
            nonlocal branch_num
            current_branch = Branch(branch_num)
            fragment_map.branches[branch_num] = current_branch
            previous_atom_data.daughter_branches.append(((bond_code, current_atom.symbol), current_branch))
            # add current branch and the bonding info to the atom that spawned current branch
            # in tuple format: ((bond info, atom symbol), branch object)
            branch_num += 1
            dfs(current_atom, previous_atom, current_branch)
            # start dfs on new branch and add info to fragment_map

        # TODO add ring detection eventually

        dfs(base_atom, None, None)
        # start a branch at base atom

        return fragment_map
    current_fragment_map = map_fragment(current_fragment, base_atom)
    print(current_fragment_map.base_atom.symbol)

    base_atom_bonds = [abbr_bond(bond) for bond in base_atom.bonded_to]
    # find base atom (i.e. highest priority atom from find_atom) from fragment
    count_base_atom_bonds = Counter(base_atom_bonds)
    # counts number of unique bonds types and how many of each base atom has

    def find_starter_atoms(atom, count_atom_bonds, atom_list):
        # searches through all atoms in molecules in total_molecule to see if they match the fragment base atom
        if atom.symbol == base_atom.symbol:
            # check to see if atom is the same element
            bond_list = [abbr_bond(bond) for bond in atom.bonded_to]
            count_bond_list = Counter(bond_list)
            for key in count_atom_bonds:
                if key not in count_bond_list or count_atom_bonds[key] > count_bond_list[key]:
                    # check 1: are there bonds types in fragment base atom that current atom doesn't have
                    # check 2: does current atom have >= the amount of each bond type compared to fragment base atom
                    # i.e. are the bonds in fragment base atom a subset of the bonds of current atom
                    return
            atom_list.append(atom)
            # if all checks passed, atom is a potential base atom and is  stored in a list

    starting_atoms = []
    # keeping track of atoms that match fragment base atom
    for atom in total_molecule.atom_list:
        find_starter_atoms(atom, count_base_atom_bonds, starting_atoms)
    print("starting atoms:")
    for atom in starting_atoms:
        print(atom.symbol)

    def check_starting_atom(starting_atom, molecule, fragment_map):

        def try_branch(branch_map, molecule_atom, previous_atom):
            atom_bonds = [abbr_bond(bond) for bond in molecule_atom.bonded_to if bond.atom != previous_atom]
            print(atom_bonds)
            print("branch_map")
            for a in branch_map.sequence:
                print(a.symbol)
                print(a.bond)
            # we are now at the first atom of the branch (which has already been confirmed to be the right element
            # need to check next bond is the same as the one in atom_bonds
            # if atom_bonds has multiple possible paths, need to apss partial branch_map to each one

        fragment_map.base_atom.daughter_branches.sort(key=lambda x: len(x[1].sequence), reverse=True)
        # this makes longer branches go first -> have to search the longest branch first
        # otherwise a shorter branch might be identified in what is actually the long branch
        # i.e. if atom has ethyl and propyl group, you could find the ethyl group where the propyl group is
        for branch in fragment_map.base_atom.daughter_branches:
            print("branch")
            for bond in starting_atom.bonded_to:
                if branch[0] == abbr_bond(bond):
                    print("potential path found from starting atom")
                    try_branch(branch[1], bond.atom, starting_atom)
            # branch[0] is a tuple containing the bond_info to the branch and the atom element of the first atom in the branch
            # looks to see if that info matches any of the bonds of the starting atom
            # if it finds a match, it does try_branch() to see if that bond path contains the branch
            # branch [1] is


    check_starting_atom(starting_atoms[0], total_molecule, current_fragment_map)

fragmentize("BrCCSi(CCCl)(CCS)CCI", "BrCCSiCCCl")

