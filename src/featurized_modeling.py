
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
df = df.loc[:, (df != 0).any(axis=0)]
df['normalized MALDI'] = df['MALDI Product Intensity']/df['MALDI Internal Standard Intensity']
df = df[df['normalized MALDI'] < 4]
df = shuffle(df, random_state=42)

df_Cu = df[df['Cu_cat'] == 1].reset_index(drop=True)
df_Ir = df[df['Ir_cat'] == 1].reset_index(drop=True)
df_Pd = df[df['Pd_cat'] == 1].reset_index(drop=True)
df_Ru = df[df['Ru_cat'] == 1].reset_index(drop=True)

data = [df_Cu, df_Ir, df_Pd, df_Ru]
catalysts = ['Cu', 'Ir', 'Pd', 'Ru']


# for i in range(len(data)):
#     df_cat = data[i]
#     X = df_cat.copy()
#     y = X.pop('EIC(+)[M+H] Product Area')
#
#     # X.drop(['MALDI Product Intensity', 'MALDI Internal Standard Intensity', 'normalized MALDI'], axis=1, inplace=True)
#     # X = df_cat.loc[:, ['MALDI Product Intensity', 'MALDI Internal Standard Intensity']]
#     kf = KFold(n_splits=5, shuffle=False)
#     fold_r_sqrs_RF = []
#     # RF_feature_importances = np.zeros(X.shape[1])
#     for train, test in kf.split(X):
#         # random forest
#         RF = RandomForestRegressor(n_estimators=500, max_depth=50, min_samples_split=3, min_samples_leaf=2, max_features=0.6, n_jobs=-1)
#         RF.fit(X.values[train], y.values[train])
#         y_pred = RF.predict(X.values[test])
#         r_sqr = r2_score(y.values[test], y_pred)
#         fold_r_sqrs_RF.append(r_sqr)
#         df_cat.loc[test, 'preds'] = pd.Series(y_pred, index=test)
#
#         # RF_feature_importances += RF.feature_importances_
# #
#         # xgboost
#         # GB = xgb.XGBRegressor(objective='reg:squarederror', learning_rate=0.01,
#         #                        max_depth=6, n_estimators=1000, subsample=0.7,
#         #                        tree_method='gpu_hist', gpu_id=0)
#         # GB.fit(X.values[train], y.values[train])
#         # y_pred = GB.predict(X.values[test])
#         # r_sqr = r2_score(y.values[test], y_pred)
#         # fold_r_sqrs_GB.append(r_sqr)
#     # avg_r_sqr_GB = np.average(np.array(fold_r_sqrs_GB))
#     # df_feature_importances = pd.DataFrame(np.array([RF_feature_importances/5]), columns=X.columns)
#     # df_feature_importances.to_csv(f'../data/feature_importances_{catalysts[i]}.csv', index=False)
#
#     avg_r_sqr_RF = np.average(np.array(fold_r_sqrs_RF))
#     print(f'catalyst: {catalysts[i]} r^2 RF: {np.round(avg_r_sqr_RF, 3)}')


#
# print(df_Cu)
# df_Cu.to_csv('../data/df_CU.csv', index=False)









# with parallel_backend('multiprocessing'):
#     for i in range(len(data)):
#         df_cat = data[i]
#         y = df_cat.pop('EIC(+)[M+H] Product Area')
#         X = df_cat.copy()
#         grid = {
#                 'n_estimators': [400, 600, 800, 1000],
#                 'max_depth': [20, 40, 60, 80, 100, None],
#                 'min_samples_split': [2, 3, 4, 5],
#                 'min_samples_leaf': [1, 2, 3],
#                 'max_features': [0.4, 0.5, 0.6, 0.7],
#                 }
#
#         model = RandomForestRegressor()
#         RF_gridsearch = GridSearchCV(estimator=model, param_grid=grid,
#                                      cv=5, verbose=0, n_jobs=-1)
#         RF_gridsearch.fit(X, y)
#         best_model = RF_gridsearch.best_estimator_
#         best_score = RF_gridsearch.best_score_
#         best_params = RF_gridsearch.best_params_
#         print(f'{catalysts[i]}, score: {best_score}, params: {best_params}')

import time
t0 = time.time()
for i in range(len(data)):
    df_cat = data[i]
    y = df_cat.pop('EIC(+)[M+H] Product Area')
    X = df_cat.copy()
    grid = {
            'learning_rate': [0.01],
            'n_estimators': [1000],
            'subsample': [0.3, 0.4, 0.5, 0.6],  # 0 to 1
            'max_depth': [3, 6, 8],  # 0 to infinite
            # 'colsample_bytree': [0.8, 0.9, 1],
            # 'colsample_bylevel': [0.8, 0.9, 1],
            # 'colsample_bynode': [0.8, 0.9, 1],
            'gamma': [0], # 0 to infinite (makes it more conservative)
            'lambda': [1], # l2 (more conservative)
            'alpha': [0], # l1 (more conservative)
            }

    model = xgb.XGBRegressor(objective='reg:squarederror', tree_method='gpu_hist', gpu_id=0)
    RF_gridsearch = GridSearchCV(estimator=model, param_grid=grid,
                                 cv=5, verbose=0)
    RF_gridsearch.fit(X, y)
    best_model = RF_gridsearch.best_estimator_
    best_score = RF_gridsearch.best_score_
    best_params = RF_gridsearch.best_params_
    print(f'{catalysts[i]}, score: {best_score}, params: {best_params}')

t1 = time.time()
print(t1-t0)