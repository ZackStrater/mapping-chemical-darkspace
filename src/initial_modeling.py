

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# Modeling only from feature available in the original dataset

df = pd.read_csv('../data/merck_data_combined.csv')

df.drop(['TWC Product Area', 'TWC Product Area%', 'MS TIC(+) Product Area%', 'MS TIC x TWC',
        'Location', 'Sample Name', 'bromide or amine', 'CHEMISTRY', 'Canonical_Smiles'], axis=1, inplace=True)
df['3/4-pyridyl-type 6 MR heterocycles (Br only)'].fillna(0, inplace=True)


df_Cu = df[df['Cu_cat'] == 1].reset_index()
df_Ir = df[df['Ir_cat'] == 1].reset_index()
df_Pd = df[df['Pd_cat'] == 1].reset_index()
df_Ru = df[df['Ru_cat'] == 1].reset_index()
data = [df_Cu, df_Ir, df_Pd, df_Ru]
catalysts = ['Cu', 'Ir', 'Pd', 'Ru']

for i in range(len(data)):
    df_cat = data[i]
    y = df_cat.pop('EIC(+)[M+H] Product Area')
    X = df_cat.copy()
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    fold_r_sqrs_RF = []
    fold_r_sqrs_GB = []
    for train, test in kf.split(X):
        # random forest
        RF = RandomForestRegressor(n_estimators=1000, max_depth=None, min_samples_split=3, min_samples_leaf=5, max_features=None)
        RF.fit(X.values[train], y.values[train])
        y_pred = RF.predict(X.values[test])
        r_sqr = r2_score(y.values[test], y_pred)
        fold_r_sqrs_RF.append(r_sqr)
        print(r_sqr)
    avg_r_sqr_RF = np.average(np.array(fold_r_sqrs_RF))
    print(f'catalyst: {catalysts[i]}\nr^2 RF: {np.round(avg_r_sqr_RF, 3)}')



