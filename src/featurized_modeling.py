
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



df = pd.read_csv('../data/featurized_data.csv')
# df = df.loc[:, (df != 0).any(axis=0)]
df['normalized MALDI'] = df['MALDI Product Intensity']/df['MALDI Internal Standard Intensity']
df = df[df['normalized MALDI'] < 4]
df = shuffle(df, random_state=42)

df_Cu = df[df['Cu_cat'] == 1].reset_index()
df_Ir = df[df['Ir_cat'] == 1].reset_index()
df_Pd = df[df['Pd_cat'] == 1].reset_index()
df_Ru = df[df['Ru_cat'] == 1].reset_index()

data = [df_Cu, df_Ir, df_Pd, df_Ru]
catalysts = ['Cu', 'Ir', 'Pd', 'Ru']

# for i in range(len(data)):
#     df_cat = data[i]
#     y = df_cat.pop('EIC(+)[M+H] Product Area')
#     X = df_cat.copy()
#     # X = df_cat.loc[:, ['MALDI Product Intensity', 'MALDI Internal Standard Intensity']]
#     kf = KFold(n_splits=5, shuffle=True, random_state=42)
#     fold_r_sqrs_RF = []
#     fold_r_sqrs_GB = []
#     for train, test in kf.split(X):
#         # random forest
#         RF = RandomForestRegressor(n_estimators=1000, max_depth=None, min_samples_split=3, min_samples_leaf=5, max_features=None, n_jobs=-1)
#         RF.fit(X.values[train], y.values[train])
#         y_pred = RF.predict(X.values[test])
#         r_sqr = r2_score(y.values[test], y_pred)
#         fold_r_sqrs_RF.append(r_sqr)
#
#         # xgboost
#         # GB = xgb.XGBRegressor(objective='reg:squarederror', learning_rate=0.01,
#         #                        max_depth=6, n_estimators=1000, subsample=0.7,
#         #                        tree_method='gpu_hist', gpu_id=0)
#         # GB.fit(X.values[train], y.values[train])
#         # y_pred = GB.predict(X.values[test])
#         # r_sqr = r2_score(y.values[test], y_pred)
#         # fold_r_sqrs_GB.append(r_sqr)
#         # print(r_sqr)
#     avg_r_sqr_RF = np.average(np.array(fold_r_sqrs_RF))
#     # avg_r_sqr_GB = np.average(np.array(fold_r_sqrs_GB))
#     print(f'catalyst: {catalysts[i]} r^2 RF: {np.round(avg_r_sqr_RF, 3)}')

from sklearn.utils import parallel_backend
with parallel_backend('multiprocessing'):
    # gridsearchcv code ...
    for i in range(len(data)):
        df_cat = data[i]
        y = df_cat.pop('EIC(+)[M+H] Product Area')
        X = df_cat.copy()
        grid = {
                'n_estimators': [750],
                'max_depth': [50],
                'min_samples_split': [3],
                'min_samples_leaf': [2],
                'max_features': [0.6],
                }

        model = RandomForestRegressor()
        RF_gridsearch = GridSearchCV(estimator=model, param_grid=grid,
                                     cv=5, verbose=0, n_jobs=-1)
        RF_gridsearch.fit(X, y)
        best_model = RF_gridsearch.best_estimator_
        best_score = RF_gridsearch.best_score_
        best_params = RF_gridsearch.best_params_
        print(f'{catalysts[i]}, score: {best_score}, params: {best_params}')

