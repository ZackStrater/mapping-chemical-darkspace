

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from find_fragments import fragmentize
from fragments_library import special_case, biomolecules, peptide_amino_acids, heterocycles, arenes, functional_groups, hydrocarbons


df = pd.read_csv('../data/merck_data_combined.csv')

y = df.pop('EIC(+)[M+H] Product Area')
X = df.copy()

X = X[0:384]
y = y[0:384]

X['3/4-pyridyl-type 6 MR heterocycles (Br only)'].fillna(0, inplace=True)
smiles_data = X['Canonical_Smiles']



for molecule in smiles_data:
    print(fragmentize(molecule, special_case, biomolecules, peptide_amino_acids, heterocycles, arenes, functional_groups, hydrocarbons))