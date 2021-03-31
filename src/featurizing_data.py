import pandas as pd
from fragments_library import special_cases, common_aromatic_heterocycles, generalized_heterocycles, secondary_amines, arenes, functional_groups, hydrocarbons, aromatic_fragments
from find_fragments import fragmentize
from smiles_to_structure import convert_to_structure, MoleculeStructure
import numpy as np
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

df = pd.read_csv('../data/merck_data_combined.csv')

df = df.loc[:, ['bromide or amine', 'Canonical_Smiles', 'MALDI Product Intensity', 'MALDI Internal Standard Intensity', 'EIC(+)[M+H] Product Area', 'Cu_cat', 'Ir_cat', 'Pd_cat', 'Ru_cat',
                'TPSA', 'cLogP', '# H-bond Acceptors', '# H-bond Donors']] # add hinderance other important columns
bromide_mask = df['bromide or amine'] == 'bromide'
amine_mask = df['bromide or amine'] == 'amine'

simple_amine_string = 'C1CNCCC1c2ccccc2'
simple_bromide_string = 'c1ncc(Br)cc1c2ccccc2'

df.loc[bromide_mask, 'amine_string'] = simple_amine_string
df.loc[bromide_mask, 'bromide_string'] = df.loc[bromide_mask, 'Canonical_Smiles']
df.loc[amine_mask, 'bromide_string'] = simple_bromide_string
df.loc[amine_mask, 'amine_string'] = df.loc[amine_mask, 'Canonical_Smiles']

libraries = [common_aromatic_heterocycles, generalized_heterocycles,secondary_amines, arenes, functional_groups, hydrocarbons, aromatic_fragments, special_cases]
fragment_feature_column_names = [frag for lib in libraries if lib != generalized_heterocycles for frag in lib]
generalized_names = ["0-5M-het", "1-5M-het", "2-5M-het", "3-5M-het", "4-5M-het",
                    "0-6M-het", "1-6M-het", "2-6M-het", "3-6M-het", "4-6M-het",
                    "0-6+5M-het", "1-6+5M-het", "2-6+5M-het", "3-6+5M-het", "4-6+5M-het", "5-6+5M-het", "6-6+5M-het",
                    "0-6+6M-het", "1-6+6M-het", "2-6+6M-het", "3-6+6M-het", "4-6+6M-het", "5-6+6M-het", "6-6+6M-het"]

fragment_feature_column_names.extend(generalized_names)
amine_feature_column_names = ['a_' + name for name in fragment_feature_column_names]
bromide_feature_column_names = ['b_' + name for name in fragment_feature_column_names]
metadata_column_names = ['mw', 'num_atoms', 'num_nitrogens', 'num_aromatic_nitrogens', 'num_non_aromatic_nitrogens', 'num_oxygens',
    'num_aromatic_oxygens', 'num_non_aromatic_oxygens', 'num_sulfurs', 'num_heteratoms', 'num_carbons', 'num_aromatic_carbons',
    'num_non_aromatic_carbons', 'num_FG', 'heterocycles', 'num_nitrogens_per_mw', 'num_oxygens_per_mw', 'num_sulfurs_per_mw',
    'num_heteratoms_per_mw', 'num_carbons_per_mw', 'num_nitrogens_per_atom', 'num_oxygens_per_atom', 'num_sulfurs_per_atom',
    'num_heteratoms_per_atom', 'num_carbons_per_atom']
amine_metadata_column_names = ['a_' + name for name in metadata_column_names]
bromide_metadata_column_names = ['b_' + name for name in metadata_column_names]

def featurize_molecules(smiles_string):
    return pd.Series(fragmentize(smiles_string, *libraries, numeric=True)[0])

def add_molecule_metadata(smiles_string):
    molecule = fragmentize(smiles_string, *libraries)[1]
    mw = molecule.mw()
    num_atoms = len(molecule.atom_list)
    num_nitrogens = 0
    num_aromatic_nitrogens = 0
    num_non_aromatic_nitrogens = 0
    num_oxygens = 0
    num_aromatic_oxygens = 0
    num_non_aromatic_oxygens = 0
    num_sulfurs = 0
    num_heteratoms = 0
    num_carbons = 0
    num_aromatic_carbons = 0
    num_non_aromatic_carbons = 0
    num_FG = len(molecule.fragments_list)
    heterocycles = 0

    heterocycle_names = ["1-5M-het", "2-5M-het", "3-5M-het", "4-5M-het",
                    "1-6M-het", "2-6M-het", "3-6M-het", "4-6M-het",
                    "1-6+5M-het", "2-6+5M-het", "3-6+5M-het", "4-6+5M-het", "5-6+5M-het", "6-6+5M-het",
                    "1-6+6M-het", "2-6+6M-het", "3-6+6M-het", "4-6+6M-het", "5-6+6M-het", "6-6+6M-het"]
    for fragment in molecule.fragments_list:
        if fragment.name in heterocycle_names or fragment.name in common_aromatic_heterocycles:
            heterocycles +=1


    for atom in molecule.atom_list:
        if atom.symbol == 'N':
            num_nitrogens += 1
            num_non_aromatic_nitrogens += 1
        if atom.symbol == 'n':
            num_nitrogens += 1
            num_aromatic_nitrogens += 1
        if atom.symbol == 'O':
            num_oxygens += 1
            num_non_aromatic_oxygens += 1
        if atom.symbol == 'o':
            num_oxygens += 1
            num_aromatic_oxygens += 1
        if atom.symbol == 'S' or atom.symbol == 's':
            num_sulfurs += 1
        if atom.heteroatom:
            num_heteratoms += 1
        if atom.symbol == 'C':
            num_carbons += 1
            num_non_aromatic_carbons += 1
        if atom.symbol == 'c':
            num_carbons += 1
            num_aromatic_carbons += 1

    num_nitrogens_per_mw = np.round(num_nitrogens/mw, 4)
    num_oxygens_per_mw = np.round(num_oxygens/mw, 4)
    num_sulfurs_per_mw = np.round(num_sulfurs/mw, 4)
    num_heteratoms_per_mw = np.round(num_heteratoms/mw, 4)
    num_carbons_per_mw = np.round(num_carbons/mw, 4)
    num_nitrogens_per_atom = np.round(num_nitrogens/num_atoms, 4)
    num_oxygens_per_atom = np.round(num_oxygens/num_atoms, 4)
    num_sulfurs_per_atom = np.round(num_sulfurs/num_atoms, 4)
    num_heteratoms_per_atom = np.round(num_heteratoms/num_atoms, 4)
    num_carbons_per_atom = np.round(num_carbons/num_atoms, 4)

    output = [mw, num_atoms, num_nitrogens, num_aromatic_nitrogens, num_non_aromatic_nitrogens, num_oxygens,
    num_aromatic_oxygens, num_non_aromatic_oxygens, num_sulfurs, num_heteratoms, num_carbons, num_aromatic_carbons,
    num_non_aromatic_carbons, num_FG, heterocycles, num_nitrogens_per_mw, num_oxygens_per_mw, num_sulfurs_per_mw,
    num_heteratoms_per_mw, num_carbons_per_mw, num_nitrogens_per_atom, num_oxygens_per_atom, num_sulfurs_per_atom,
    num_heteratoms_per_atom, num_carbons_per_atom]

    return pd.Series(output)

def get_coupling_fragments(amine_smiles_string, bromide_smiles_string):
    amine = fragmentize(amine_smiles_string, *libraries)[1]
    bromide = fragmentize(amine_smiles_string, *libraries)[1]
    amine_fragment_names = ["piperazine", "morpholine", "pyrrolidine", "pyrazolidine", "piperidine", "amine(aryl)", "primary amine", "secondary amine"]
    bromide_fragment_names = []
    # get het and get amine
    # add arylamines mb
    # get 2d adjacency matrix





print(df)
df[amine_feature_column_names] = df['amine_string'].apply(lambda x: featurize_molecules(x))
df[amine_metadata_column_names] = df['amine_string'].apply(lambda x: add_molecule_metadata(x))
df[bromide_feature_column_names] = df['bromide_string'].apply(lambda x: featurize_molecules(x))
df[bromide_metadata_column_names] = df['bromide_string'].apply(lambda x: add_molecule_metadata(x))

df.drop(['bromide or amine', 'Canonical_Smiles', 'amine_string', 'bromide_string'], axis=1, inplace=True)
print(df)
df.to_csv('../data/featurized_data.csv', index=False)