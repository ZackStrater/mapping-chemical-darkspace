import pandas as pd
from fragments_library import special_cases, common_aromatic_heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments
from find_fragments import fragmentize
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

df = pd.read_csv('../data/merck_data_combined.csv')

y = df.pop('EIC(+)[M+H] Product Area')
X = df.loc[:, ['Canonical_Smiles', 'MALDI Product Intensity', 'MALDI Internal Standard Intensity', 'Cu_cat', 'Ir_cat', 'Pd_cat', 'Ru_cat']]

libraries = [common_aromatic_heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments, special_cases]
fragment_column_names = [frag for lib in libraries for frag in lib]
def featurize_molecules(smiles_string):
    return pd.Series(fragmentize(smiles_string, *libraries, numeric=True)[0])


generalized_names = ["0-5M-het", "1-5M-het", "2-5M-het", "3-5M-het", "4-5M-het",
                    "0-6M-het", "1-6M-het", "2-6M-het", "3-6M-het", "4-6M-het",
                    "0-6+5M-het", "1-6+5M-het", "2-6+5M-het", "3-6+5M-het", "4-6+5M-het", "5-6+5M-het", "6-6+5M-het",
                    "0-6+6M-het", "1-6+6M-het", "2-6+6M-het", "3-6+6M-het", "4-6+6M-het", "5-6+6M-het", "6-6+6M-het"]

print(fragment_column_names)
print(len(fragment_column_names))
print((fragmentize('C', *libraries, numeric=True)[0]))
X[fragment_column_names] = X['Canonical_Smiles'].apply(lambda x: featurize_molecules(x))
print(X)