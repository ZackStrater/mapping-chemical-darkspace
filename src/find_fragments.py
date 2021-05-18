
from src.smiles_to_structure import convert_to_structure, MoleculeStructure, Fragment
from collections import Counter
from termcolor import cprint
from src.fragments_library import special_cases, biomolecules, peptide_amino_acids, heterocycles, \
    common_aromatic_heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments


'''The find_fragment method allows a substructure search of a given chemical fragment within a molecule 
(subgraph isomorphism).  Takes 
  
The fragmentize method allows a library of fragments to be searched for within a molecule via the find_fragment method.  
Takes the SMILES string of a molecule and fragment libraries as *args as inputs.

The fragment library should be ordered hierarchically, with more complex fragments being searched first.  Atoms found in 
any given substructure search are marked as discovered and are not used in further substructure searches (unless 
specified as "phantom atoms", see fragments_library).  
The method will return a list of names of the fragments found and the labeled molecular structure as a tuple 
-> (fragment_names_list, labeled_structure)

If numeric=True, will return a vector with the count of the number of each fragment found and the labeled molecular 
structure as a tuple
-> (fragment_vector, labeled_structure)
'''


class AtomData:
    def __init__(self, symbol):
        self.symbol = symbol
        self.bond = None
        self.daughter_branches = []
        self.ring_closures = set()
        self.phantom_bonds = None
        self.phantom_atom = False


class Branch:
    def __init__(self, bond):
        self.bond = bond
        self.sequence = []


def abbr_bond(bond):
    return bond.bond_code, bond.atom.symbol


# bond_atom_info is the atom_info of the atom for the bond being checked (to see if it's a phantom bond/discovered)
def check_bond(bond, map_bond, bond_atom_info):

    # if the bond.atom is discovered already, it should give back false (can't retread over discovered atoms)
    # unless the atom_info for that atom shows that the bond.atom should be a phantom atom (which can be searched for in discovered atoms)
    if bond.atom.discovered:
        if not bond_atom_info.phantom_atom:
            return False

    if abbr_bond(bond) == map_bond:
        return True
    # need to cover (correct, R) (correct, Q) (9, R) (9, Q) (9, correct)
    elif bond.bond_code == map_bond[0] or map_bond[0] == 9: # covers (correct, R) (correct, Q) (9, R) (9, Q)
        if map_bond[1] == "R":
            return True
        elif map_bond[1] == "Q" and bond.atom.heteroatom:
            return True
    elif bond.atom.symbol == map_bond[1]:
        if map_bond[0] == 9:
            return True
    else:
        return False


def calc_branch_length(branch):
    branch_length = 0

    def add_daughter_branch_length(daughter_branch):
        nonlocal branch_length
        branch_length += len(daughter_branch.sequence)
        if len(daughter_branch.sequence[-1].daughter_branches) > 0:
            for b in daughter_branch.sequence[-1].daughter_branches:
                add_daughter_branch_length(b)

    add_daughter_branch_length(branch)
    return branch_length


def find_fragment(fragment_string, molecule_string, fragment_name, structure=None, verbose=False):
    verbose_bin = []

    if structure:
        molecule_structure = structure
    else:
        molecule_structure = convert_to_structure(MoleculeStructure(), molecule_string)
    fragment_structure = convert_to_structure(MoleculeStructure(), fragment_string)

    def find_anchor_atom(fragment):
        for ele in ["Si", "P", "p", "S", "s", "I", "Br", "Cl", "F", "B", "b", "O", "o", "N", "n", "C", "c", "R"]:
            for atom in fragment.atom_list:
                if atom.symbol == ele:
                    return atom

    fragment_anchor_atom = find_anchor_atom(fragment_structure)
    # the actual atom object of highest priority in the fragment structure

    def is_potential_anchor(atom, fragment_anchor_atom, atom_list):
        # searches through all atoms in molecules in total_molecule to see if they match the fragment base atom
        # atom -> current atom its checking
        # atom_list is list where potential anchor atoms are stored
        # fragment_anchor_atom is the actual atom object from the fragment structure

        if atom.discovered and not fragment_anchor_atom.phantom_atom:  # if fragment_anchor atom is a phantom atom, it can use discovered atoms as potential anchors
            return
        # atom has already been used to find a fragment

        if atom.symbol != fragment_anchor_atom.symbol and fragment_anchor_atom.symbol != 'R':  # TODO what about if anchor atom is Q!??
            return
            # check to see if atom is the same element

        fragment_anchor_atom_bonds = Counter([abbr_bond(bond) for bond in fragment_anchor_atom.bonded_to])
        # count bonds from anchor atom

        atom_bonds = Counter([abbr_bond(bond) for bond in atom.bonded_to])

        # count bonds in potential anchor atom where the bond's atom haven't been discovered yet (as we won't be able to use those bonds)
        for key in fragment_anchor_atom_bonds:
            if key[1] != "R" and key[1] != "Q" and key[0] != 9:  # TODO better way to do this???
                if key not in atom_bonds or fragment_anchor_atom_bonds[key] > atom_bonds[key]:
                    # check 1: are there bonds types in fragment base atom that current atom doesn't have
                    # check 2: does current atom have >= the amount of each bond type compared to fragment base atom
                    # i.e. are the bonds in fragment anchor atom a subset of the bonds of current atom
                    return
        atom_list.append(atom)
        # if all checks passed, atom is a potential base atom and is  stored in a list

    potential_anchor_atoms = []
    # keeping track of atoms that match fragment base atom
    for atom in molecule_structure.atom_list:
        is_potential_anchor(atom, fragment_anchor_atom, potential_anchor_atoms)

    if potential_anchor_atoms == []:
        verbose_bin.append("no anchor atoms found")
        return 0
    else:
       verbose_bin.append("potential anchor atoms: ")
       for atom in potential_anchor_atoms:
           verbose_bin.append(atom.symbol)
           for bond in atom.bonded_to:
               verbose_bin.append(abbr_bond(bond))

    def map_fragment(fragment, anchor_atom):

        visited = {}
        for atom in fragment.atom_list:
            visited[atom] = False
        # keeps track of which atoms have been visited

        atom_info_dict = {}
        # links the molecule_atom and the atom_info representing that atom, used to pass ring_closure info to map

        ring_closure_counter = 1

        def traverse(current_atom, previous_atom, current_branch):
            visited[current_atom] = True

            current_atom_data = AtomData(current_atom.symbol)
            # data object for current atom

            # atom_data will reflect that the atom is a phantom_atom
            if current_atom.phantom_atom:
                current_atom_data.phantom_atom = True

            atom_info_dict[current_atom] = current_atom_data

            if current_branch:
                current_branch.sequence.append(current_atom_data)
                # append atom info to branch sequence
                # if current_branch b/c first atom does not have a branch

            current_atom_data.phantom_bonds = current_atom.phantom_bonds

            unchecked_bonds = [bond for bond in current_atom.bonded_to if bond.atom != previous_atom]

            nonlocal ring_closure_counter

            # if more than 1 unchecked bonds (i.e. a branch point), create new branch for each unchecked bond
            if len(unchecked_bonds) > 1:
                for bond in unchecked_bonds:
                    if not visited[bond.atom]:
                        verbose_bin.append("new branch")
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)
                    elif not bool(current_atom_data.ring_closures & atom_info_dict[bond.atom].ring_closures):
                        # if visited[bond.atom], we are at a ring closure
                        # this bool sees if the atom_info of these two atoms (current atom and the atom its bonded to) share any values (& operator)
                        # if they do, this ring closure has already been documented and we don't want to double count it
                        verbose_bin.append("ring closure")
                        current_atom_data.ring_closures.add((ring_closure_counter, bond.bond_code))
                        atom_info_dict[bond.atom].ring_closures.add((ring_closure_counter, bond.bond_code))
                        # add matching values to each atom_info.ring_closure
                        # ring closure data in format (ring closure #, bond_code)
                        ring_closure_counter += 1

            # if a contiguous section of branch, add bond info
            elif len(unchecked_bonds) == 1:
                if current_branch:
                    if not visited[unchecked_bonds[0].atom]:
                        verbose_bin.append("continue branch")
                        current_atom_data.bond = abbr_bond(unchecked_bonds[0])
                        traverse(unchecked_bonds[0].atom, current_atom, current_branch)
                    elif not bool(current_atom_data.ring_closures & atom_info_dict[unchecked_bonds[0].atom].ring_closures):
                        verbose_bin.append("ring closure")
                        current_atom_data.ring_closures.add((ring_closure_counter, unchecked_bonds[0].bond_code))
                        atom_info_dict[unchecked_bonds[0].atom].ring_closures.add((ring_closure_counter, unchecked_bonds[0].bond_code))
                        ring_closure_counter += 1
                        # same as above
                else:
                    verbose_bin.append("new branch")
                    for bond in unchecked_bonds:
                        new_branch(bond.atom, current_atom, current_atom_data, bond.bond_code)
                    # if the anchor atom only has 1 bond, need to start a branch

            else:
                verbose_bin.append("end point")

            if not current_branch:
                return current_atom_data
                # this returns anchor atom to the map_fragment function as the anchor atom is not spawned from a branch @

        def new_branch(current_atom, previous_atom, previous_atom_data, bond_code):
            current_branch = Branch((bond_code, current_atom.symbol))
            # create new branch with bonding info to first atom in branch
            previous_atom_data.daughter_branches.append(current_branch)
            # add new branch to the atom which spawned it
            traverse(current_atom, previous_atom, current_branch)
            # start traverse on first atom in branch
            # need to pass previous_atom in order to not travel backwards

        return traverse(anchor_atom, None, None)
        # starts process of mapping fragment, but also returns the anchor atom

    anchored_fragment_map = map_fragment(fragment_structure, fragment_anchor_atom)
    # the map base is the atom_data representation of the anchor atom
    # the rest of the map is stored in the daughter branches

    def expand_map(anchor_atom):
        verbose_bin.append("anchor atom")
        verbose_bin.append(anchor_atom.symbol)
        if len(anchor_atom.ring_closures) > 0:
            verbose_bin.append("ring closures:")
            for num in anchor_atom.ring_closures:
                verbose_bin.append(num)
        if anchor_atom.phantom_bonds:
            verbose_bin.append(f"phantom bonds = {anchor_atom.phantom_bonds}")

        def expand_branch_point(atom_map):
            for branch in atom_map.daughter_branches:
                verbose_bin.append("branch:")
                verbose_bin.append(f"branch length: {len(branch.sequence)}")
                verbose_bin.append(f"total branch length: {calc_branch_length(branch)}")
                verbose_bin.append(f"bond to branch: {branch.bond}")
                for atom_info in branch.sequence:
                    verbose_bin.append(atom_info.symbol)
                    if len(atom_info.ring_closures) > 0:
                        verbose_bin.append("ring closures:")
                        for num in atom_info.ring_closures:
                            verbose_bin.append(num)
                    if atom_info.phantom_bonds:
                        verbose_bin.append(f"phantom bonds = {atom_info.phantom_bonds}")
                    if atom_info.bond:
                        verbose_bin.append(atom_info.bond)
                    if len(atom_info.daughter_branches) > 0:
                        verbose_bin.append("branch point")
                        expand_branch_point(atom_info)
        expand_branch_point(anchor_atom)


    verbose_bin.append("\nexpanded map:\n")
    expand_map(anchored_fragment_map)

    def check_anchor_atom(potential_anchor_atom, fragment_map):
        molecule_atoms = {potential_anchor_atom}
        # list to keep track of which atoms in the molecule constitute a matched fragment

        currently_visited = {potential_anchor_atom: fragment_map}
        # dictionary that keeps track of which atoms have been used to find the fragment at any given step

        def check_branch_point(current_molecule_atom, previous_molecule_atom, map_atom_info, branch_atoms):
            if map_atom_info.phantom_bonds:
                bond_num = len(current_molecule_atom.bonded_to)
                if bond_num != map_atom_info.phantom_bonds:
                    verbose_bin.append("wrong amount of phantom bonds")
                    return False
            # phantom_bonds is a way to ensure the current atom is bonded to the specified number of atoms
            # note that phantom bonds includes any bonds for the current molecule_atom, including those to atoms that are "discovered"

            branch_point_atoms = set()
            nonlocal currently_visited
            verbose_bin.append("I'm trying a branch point")
            map_atom_info.daughter_branches.sort(key=calc_branch_length, reverse=True)
            # this makes longer branches go first -> have to search the longest branch first
            # otherwise a shorter branch might be identified in what is actually the long branch
            # i.e. if atom has ethyl and propyl group, you could find the ethyl group where the propyl group is and then be unable to find propyl group
            # also important - need to caculate the total branch length (including length of all its daughter branches)

            verbose_bin.append("branch point bonds check")
            unchecked_bonds = Counter([abbr_bond(bond) for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom])
            fragment_branch_point_bonds = Counter([branch.bond for branch in map_atom_info.daughter_branches])
            verbose_bin.append(unchecked_bonds)
            verbose_bin.append(fragment_branch_point_bonds)
            # subset check on branch point, just to make sure current atom has all the bonds the fragment branchpoint has
            for key in fragment_branch_point_bonds:
                if key[1] != "R" and key[1] != "Q" and key[0] != 9:  # TODO better way to do this?
                    if key not in unchecked_bonds or fragment_branch_point_bonds[key] > unchecked_bonds[key]:
                        verbose_bin.append("branch point doesn't contain necessary bonds")
                        return False

            branch_check = {}
            for branch in map_atom_info.daughter_branches:
                branch_check[branch] = False
            # set all branches to unconfirmed

            trial_paths = [bond for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom]
            # available routes to see if branch matches
            checked_paths = []
            # keeps track of bonds that have been used to successfully identify branches
            for branch in map_atom_info.daughter_branches:
                # take each branch
                for bond in trial_paths:
                    if check_bond(bond, branch.bond, branch.sequence[0]) and bond not in checked_paths and bond.atom not in currently_visited:
                        # if the bond to the branch matches the current bond (and the current bond hasn't already been used to identify a branch):
                        if try_branch(branch.sequence, bond.atom, current_molecule_atom, branch_point_atoms):
                            # test to see if the current branch works on this bond path
                            verbose_bin.append("branch successful")
                            branch_check[branch] = True
                            checked_paths.append(bond)
                            # if true, the branch was successfully found, turn branch to True in branch check
                            # add bond used to checked_paths so it isn't used for further branches
                            break
                            # need to stop the for loop so it doesn't search the matched branch in further trial_paths
                        else:
                            verbose_bin.append("branch not successful")

            if all(value is True for value in branch_check.values()):
                verbose_bin.append("branch point match")
                if branch_atoms:
                    branch_atoms.update(branch_point_atoms)
                else:
                    molecule_atoms.update(branch_point_atoms)
                    # first branch point does not have a branch that spawned it
                return True
                # if all branches have been found, they will be True in branch_check, branch point is a match, return True
            else:
                verbose_bin.append("branch point not matched")
                return False
                # one or more branches were found, branch point wasn't a match, return False

        def try_branch(branch_sequence, current_molecule_atom, previous_molecule_atom, branch_point_atoms):
            branch_atoms = set()
            verbose_bin.append("I'm trying a branch!")
            if check_atom_bonds(current_molecule_atom, previous_molecule_atom, branch_sequence, 0, branch_atoms):
                branch_point_atoms.update(branch_atoms)
                # add branch_atoms to branch point_atoms
                return True
            else:
                nonlocal currently_visited
                for a in branch_atoms:
                    currently_visited.pop(a)

        def check_ring_closure(current_atom, atom_info):  # already check if atom_info has ring closure
            ring_closures = set()  # all the ring closure numbers in currently_visited
            for value in currently_visited.values():
                ring_closures.update(value.ring_closures)

            for closure_info in atom_info.ring_closures:  # checking each ring closure atom has
                verbose_bin.append("ring closure")
                verbose_bin.append(closure_info)
                if closure_info in ring_closures:  # we've already hit the other half of the ring closure
                    for key in currently_visited:  # looking for matching ring closure
                        if closure_info in currently_visited[key].ring_closures:  # matched ring closure, key = atom it should be bonded to
                            ring_closure_partner = key

                            if ring_closure_partner in [bond.atom for bond in current_atom.bonded_to]:
                                ring_closure_bond = None
                                for bond in current_atom.bonded_to:
                                    if bond.atom == ring_closure_partner:
                                        ring_closure_bond = bond
                                if ring_closure_bond.bond_code != closure_info[1] and closure_info[1] != 9:
                                    verbose_bin.append("closure bond incorrect")
                                    return False
                            else:
                                verbose_bin.append("atom not bonded to correct closure partner")
                                return False
                else:
                    return True
                # first time encountering that ring closure number, don't need to do any further checks
            verbose_bin.append("all ring closures acounted for")
            return True

        def check_atom_bonds(current_molecule_atom, previous_molecule_atom, branch_sequence, index, branch_atoms):
            verbose_bin.append("checking atom")
            verbose_bin.append(current_molecule_atom.symbol)
            nonlocal currently_visited
            current_atom_info = branch_sequence[index]

            if current_atom_info.phantom_bonds:
                bond_num = len(current_molecule_atom.bonded_to)
                if bond_num != current_atom_info.phantom_bonds:
                    verbose_bin.append("wrong amount of phantom bonds")
                    return False

            if current_atom_info.ring_closures:
                if not check_ring_closure(current_molecule_atom, current_atom_info):
                    return False

            currently_visited[current_molecule_atom] = current_atom_info
            for a in currently_visited:
                verbose_bin.append(a.symbol)
                if currently_visited[a].ring_closures:
                    verbose_bin.append(currently_visited[a].ring_closures)
            verbose_bin.append("\n")

            branch_atoms.add(current_molecule_atom)
            if len(current_atom_info.daughter_branches) > 0:
                # atom is branch point and need to check branches
                return check_branch_point(current_molecule_atom, previous_molecule_atom, current_atom_info, branch_atoms)
            else:
                # atom is either an endpoint or contiguous segment:
                if not current_atom_info.bond:
                    verbose_bin.append("reached branch end")

                    # if no bond data, means we have matched the entire branch, return True
                    return True
                else:
                    # else: this is a contiguous segment look for appropriate bonds
                    unchecked_bonds = [bond for bond in current_molecule_atom.bonded_to if bond.atom != previous_molecule_atom]
                    if len(unchecked_bonds) == 0:
                        verbose_bin.append("branch ended too early")
                        # actual branch has ended, but map says there should be another atom bonded here, therefore return False
                        return False
                    elif len(unchecked_bonds) == 1:
                        # actual molecule only has a contiguous segment here
                        verbose_bin.append(current_atom_info.bond)
                        verbose_bin.append(abbr_bond(unchecked_bonds[0]))
                        if check_bond(unchecked_bonds[0], current_atom_info.bond, branch_sequence[index + 1]) and unchecked_bonds[0].atom not in currently_visited:
                            return check_atom_bonds(unchecked_bonds[0].atom, current_molecule_atom, branch_sequence, index + 1, branch_atoms)  # check next atom
                            # all branches should either return a function, True, or False.  All child functions should do the same
                            # uncheck_bonds[0].atom becomes new current_molecule_atom, current_molecule_atom becomes previous_molecule_atom
                            # also pass the branch sequence and the index of the next atom_info in the branch
                        else:
                            verbose_bin.append("wrong bond or already visited")
                            return False
                            # the next atom in the branch either has the wrong bond or atom symbol
                    else:
                        # there are multiple possible paths branch could go
                        verbose_bin.append("checking multiple paths for branch")
                        # check all ways
                        for bond in unchecked_bonds:
                            if current_atom_info.bond != abbr_bond(bond):  # this is purely for seeing what's happening
                                verbose_bin.append(abbr_bond(bond))
                                verbose_bin.append(current_atom_info.bond)
                            if check_bond(bond, current_atom_info.bond, branch_sequence[index + 1]) and bond.atom not in currently_visited:
                                verbose_bin.append(abbr_bond(bond))
                                verbose_bin.append(current_atom_info.bond)
                                # looks at all possible ways that match the correct bond
                                midway_fork = set()
                                # need to separate the branch_atoms here since we don't know if any of the paths will work
                                if check_atom_bonds(bond.atom, current_molecule_atom, branch_sequence, index + 1, midway_fork):
                                    branch_atoms.update(midway_fork)
                                    # if one of the paths works, add all the atoms from the midway_fork "branch"
                                    return True
                                    # return True if any of the paths work (also returns first found)
                                else:
                                    for a in midway_fork:
                                        currently_visited.pop(a)
                                    # if midway_fork doesn't work, need to remove those atoms from currently_visited
                        return False
                        # if this else clause didn't return True (i.e. none of the paths succeeded)
                        # then none of the paths are productive, return false

        if check_branch_point(potential_anchor_atom, None, fragment_map, None):
            verbose_bin.append("phantom atom check")
            for atom in currently_visited:
                if not currently_visited[atom].phantom_atom:
                    # using currently visited to see if the atom_data that was used to find that atom was marked as a phantom_atom
                    # if it is a phantom atom, don't mark as discovered
                    atom.discovered = True
                else:
                    verbose_bin.append("this atom should not be counted")

            # running checks
            verbose_bin.append(f"number of atoms in fragment: {len(molecule_atoms)}")
            molecule_atoms2 = [atom.symbol for atom in molecule_atoms]
            molecule_atoms2phantom = [atom.phantom_atom for atom in molecule_atoms]
            if (len(molecule_atoms)) != len(currently_visited):
                verbose_bin.append("error in number of atoms found")
            for atom in molecule_atoms:
                if atom not in currently_visited:
                    verbose_bin.append("error: descrepancy between currently_visited and molecule_atoms")

            for atom in molecule_atoms:
                verbose_bin.append(atom.symbol)
            verbose_bin.append("matched fragment to anchor atom")
            return Fragment(fragment_name, list(molecule_atoms)) # TODO currently this includes atoms that were found via phantom atoms (unclear if this is wanted behavior)
        else:
            verbose_bin.append("anchor atom not matched to fragment")
            return False
        # start from check_branch point on the potential anchor atom
        # the anchor atom in map is treated as a branch point, even if it only has 1 branch



    fragment_counter = 0
    fragments_identified = []
    for atom in potential_anchor_atoms:
        verbose_bin.append("checking anchor atom")
        for bond in atom.bonded_to:
            verbose_bin.append(abbr_bond(bond))
        is_found_fragment = check_anchor_atom(atom, anchored_fragment_map)
        if is_found_fragment:
            # add atoms found to fragment
            fragment_counter += 1
            fragments_identified.append(is_found_fragment)

    verbose_bin.append(f"\nnumber of fragments found: {fragment_counter}")
    if verbose:
        for item in verbose_bin:
            print(item)

    return fragments_identified



def fragmentize(molecule_string, *fragment_libraries, numeric=False, verbose=False):

    molecule_structure = convert_to_structure(MoleculeStructure(), molecule_string)
    fragments = []
    fragment_names = []
    fragments_counter = []

    generalized_heterocycles_found = []
    for lib in fragment_libraries:
        if lib != generalized_heterocycles:
            for frag in lib:
                frag_num = 0
                for frag_res_structure in lib[frag]:
                    frag_res_found = find_fragment(frag_res_structure, None, frag, structure=molecule_structure, verbose=verbose)
                    if frag_res_found:
                        frag_num += len(frag_res_found)
                        # can find multiples of a fragment
                        for f in frag_res_found:
                            fragments.append(f)
                            fragment_names.append(f.name)
                fragments_counter.append(frag_num)
        # for generalized heterocycles
        else:
            for frag in lib:
                for frag_res_structure in lib[frag]:
                    frag_res_found = find_fragment(frag_res_structure, None, frag, structure=molecule_structure, verbose=verbose)
                    if frag_res_found:
                        for f in frag_res_found:
                            f.generalize_heterocycle_name()
                            generalized_heterocycles_found.append(f)

    # possible varieties of generalized heterocycles
    # name format: X-Y+ZM-het where X is number of heteroatoms, Y is the number of atoms in the ring and
    # Z is the number of atoms in the fused ring
    generalized_heterocycles_names = ["0-5M-het", "1-5M-het", "2-5M-het", "3-5M-het", "4-5M-het",
                                      "0-6M-het", "1-6M-het", "2-6M-het", "3-6M-het", "4-6M-het",
                                      "0-6+5M-het", "1-6+5M-het", "2-6+5M-het", "3-6+5M-het", "4-6+5M-het", "5-6+5M-het", "6-6+5M-het",
                                      "0-6+6M-het", "1-6+6M-het", "2-6+6M-het", "3-6+6M-het", "4-6+6M-het", "5-6+6M-het", "6-6+6M-het"]


    generalized_heterocycles_found_dict = {k:0 for k in generalized_heterocycles_names}
    for heterocycle in generalized_heterocycles_found:
        generalized_heterocycles_found_dict[heterocycle.name] += 1
        fragments.append(heterocycle)
        fragment_names.append(heterocycle.name)
    for key in generalized_heterocycles_found_dict:
        fragments_counter.append(generalized_heterocycles_found_dict[key])

    molecule_structure.fragments_list = fragments
    atoms_not_discovered = 0
    for atom in molecule_structure.atom_list:
        if not atom.discovered:
            atoms_not_discovered += 1
    if atoms_not_discovered > 0:
        # total_frags = 0
        # for lib in fragment_libraries:
        #     total_frags += len(lib)
        cprint(f"atoms not found: {atoms_not_discovered}", "red")
        cprint(molecule_string, 'red')
        # return ["NA" for _ in range(total_frags)]
    if numeric:
        return fragments_counter, molecule_structure
    else:
        return fragment_names, molecule_structure


