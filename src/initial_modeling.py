

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


df = pd.read_csv('../data/merck_data_combined.csv')

y = df.pop('EIC(+)[M+H] Product Area')
X = df.copy()
X.drop(['TWC Product Area', 'TWC Product Area%', 'MS TIC(+) Product Area%', 'MS TIC x TWC',
        'Location', 'Sample Name', 'bromide or amine', 'CHEMISTRY', 'Canonical_Smiles'], axis=1, inplace=True)

X_train, X_test, y_train, y_test = train_test_split(X, y, )