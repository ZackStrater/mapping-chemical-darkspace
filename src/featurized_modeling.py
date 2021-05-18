
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.model_selection import GridSearchCV
from sklearn.utils import parallel_backend, shuffle
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt


# Modeling using features derived from featurizing_data

df = pd.read_csv('../data/featurized_data.csv')
df = df.loc[:, (df != 0).any(axis=0)]
df = df[df['normalized MALDI'] < 4]
df = shuffle(df, random_state=42)

# Reactions with each catalysts are modeled separately
df_Cu = df[df['Cu_cat'] == 1].reset_index(drop=True)
df_Ir = df[df['Ir_cat'] == 1].reset_index(drop=True)
df_Pd = df[df['Pd_cat'] == 1].reset_index(drop=True)
df_Ru = df[df['Ru_cat'] == 1].reset_index(drop=True)

data = [df_Cu, df_Ir, df_Pd, df_Ru]
catalysts = ['Cu', 'Ir', 'Pd', 'Ru']


for i in range(len(data)):
    df_cat = data[i]
    X = df_cat.copy()
    X.drop(['amine_string', 'bromide_string'], axis=1, inplace=True)
    # X = df_cat.loc[:, 'MALDI Product Intensity': 'normalized MALDI']
    # X.drop(['MALDI Product Intensity', 'MALDI Internal Standard Intensity', 'normalized MALDI'], axis=1, inplace=True)
    # X = df_cat.loc[:, 'MALDI Product Intensity': 'b_num_carbons_per_atom']
    y = X.pop('EIC(+)[M+H] Product Area')

    kf = KFold(n_splits=5, shuffle=False)
    fold_r_sqrs_RF = []
    RF_feature_importances = np.zeros(X.shape[1])
    for train, test in kf.split(X):
        # random forest
        RF = RandomForestRegressor(n_estimators=500, max_depth=50, min_samples_split=3, min_samples_leaf=2, max_features=0.6, n_jobs=-1)
        RF.fit(X.values[train], y.values[train])
        y_pred = RF.predict(X.values[test])
        r_sqr = r2_score(y.values[test], y_pred)
        fold_r_sqrs_RF.append(r_sqr)
        df_cat.loc[test, 'preds'] = pd.Series(y_pred, index=test)
        RF_feature_importances += RF.feature_importances_

    avg_r_sqr_RF = np.mean(np.array(fold_r_sqrs_RF))
    print(f'catalyst: {catalysts[i]} r^2 RF: {np.round(avg_r_sqr_RF, 3)}')


df_predictions = pd.concat([df_Cu, df_Ir, df_Pd, df_Ru], axis=0)
df_predictions.to_csv('../data/predictions.csv', index=False)




