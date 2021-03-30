import pandas as pd
from fragments_library import special_cases, common_aromatic_heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments
from find_fragments import fragmentize
import numpy as np
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

df = pd.read_csv('../data/merck_data_combined.csv')

y = df.pop('EIC(+)[M+H] Product Area')
X = df.loc[:, ['bromide or amine', 'Canonical_Smiles', 'MALDI Product Intensity', 'MALDI Internal Standard Intensity', 'Cu_cat', 'Ir_cat', 'Pd_cat', 'Ru_cat']]
bromide_mask = X['bromide or amine'] == 'bromide'
amine_mask = X['bromide or amine'] == 'amine'

simple_amine_string = 'C1CNCCC1c2ccccc2'
simple_bromide_string = 'c1ncc(Br)cc1c2ccccc2'

X.loc[bromide_mask, 'amine_string'] = simple_amine_string
X.loc[bromide_mask, 'bromide_string'] = X.loc[bromide_mask, 'Canonical_Smiles']
X.loc[amine_mask, 'bromide_string'] = simple_bromide_string
X.loc[amine_mask, 'amine_string'] = X.loc[amine_mask, 'Canonical_Smiles']

libraries = [common_aromatic_heterocycles, generalized_heterocycles, arenes, functional_groups, hydrocarbons, aromatic_fragments, special_cases]
fragment_feature_column_names = [frag for lib in libraries if lib != generalized_heterocycles for frag in lib]
generalized_names = ["0-5M-het", "1-5M-het", "2-5M-het", "3-5M-het", "4-5M-het",
                    "0-6M-het", "1-6M-het", "2-6M-het", "3-6M-het", "4-6M-het",
                    "0-6+5M-het", "1-6+5M-het", "2-6+5M-het", "3-6+5M-het", "4-6+5M-het", "5-6+5M-het", "6-6+5M-het",
                    "0-6+6M-het", "1-6+6M-het", "2-6+6M-het", "3-6+6M-het", "4-6+6M-het", "5-6+6M-het", "6-6+6M-het"]

fragment_feature_column_names.extend(generalized_names)
amine_feature_column_names = ['a_' + name for name in fragment_feature_column_names]
bromide_feature_column_names = ['b_' + name for name in fragment_feature_column_names]



def featurize_molecules(smiles_string):
    return pd.Series(fragmentize(smiles_string, *libraries, numeric=True)[0])

print(X)
X[amine_feature_column_names] = X['amine_string'].apply(lambda x: featurize_molecules(x))
X[bromide_feature_column_names] = X['bromide_string'].apply(lambda x: featurize_molecules(x))
X.drop(['bromide or amine', 'Canonical_Smiles', 'amine_string', 'bromide_string'], axis=1, inplace=True)
print(X)
X.to_csv('../data/featurized_Xdata.csv', index=False)