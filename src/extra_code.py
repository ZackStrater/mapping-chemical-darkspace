

## FROM FEATURIZED MODELING to compare outputs of same reaction
# # merging out of sample predictions
# df_Cu_pred = df_Cu.loc[:, ['normalized MALDI', 'MALDI Product Intensity', 'EIC(+)[M+H] Product Area', 'preds', 'amine_string', 'bromide_string']]
# df_Cu_pred.columns = ['Cu Normalized MALDI', 'Cu MALDI Product Intensity', 'Cu EIC(+)[M+H] Product Area', 'Cu preds', 'amine_string', 'bromide_string']
# df_Ir_pred = df_Ir.loc[:, ['normalized MALDI', 'MALDI Product Intensity', 'EIC(+)[M+H] Product Area', 'preds', 'amine_string', 'bromide_string']]
# df_Ir_pred.columns = ['Ir Normalized MALDI', 'Ir MALDI Product Intensity', 'Ir EIC(+)[M+H] Product Area', 'Ir preds', 'amine_string', 'bromide_string']
# df_Pd_pred = df_Pd.loc[:, ['normalized MALDI', 'MALDI Product Intensity', 'EIC(+)[M+H] Product Area', 'preds', 'amine_string', 'bromide_string']]
# df_Pd_pred.columns = ['Pd Normalized MALDI', 'Pd MALDI Product Intensity', 'Pd EIC(+)[M+H] Product Area', 'Pd preds', 'amine_string', 'bromide_string']
# df_Ru_pred = df_Ru.loc[:, ['normalized MALDI', 'MALDI Product Intensity', 'EIC(+)[M+H] Product Area', 'preds', 'amine_string', 'bromide_string']]
# df_Ru_pred.columns = ['Ru Normalized MALDI', 'Ru MALDI Product Intensity', 'Ru EIC(+)[M+H] Product Area', 'Ru preds', 'amine_string', 'bromide_string']
# dfs = [df_Cu_pred, df_Ir_pred, df_Pd_pred, df_Ru_pred]
# from functools import reduce
#
# df_predictions = reduce(lambda left, right: pd.merge(left, right, on=['amine_string', 'bromide_string']), dfs)
#
# df_predictions.to_csv('../data/predictions_merged.csv', index=False)

# # FROM FEATURIZED MODELING to compare feature importances
# features importances while modeling
# df_feature_importances = pd.DataFrame(np.array([RF_feature_importances/5]), columns=X.columns)
# maldi_feature_importances = df_feature_importances.loc[:, 'MALDI Product Intensity':'normalized MALDI']
# meta_data_importances = df_feature_importances.loc[:, 'TPSA':'b_num_carbons_per_atom']
# total_fragment_importances = df_feature_importances.loc[:, 'a_indole':'bromide_het_neighbor_0-6M-het']

# print(f"{catalysts[i]} importances:")
# print(f"maldi: {maldi_feature_importances.sum(axis=1)}\nmetadata: {meta_data_importances.sum(axis=1)}\nfragments: {total_fragment_importances.sum(axis=1)}")
# print(f"metadata: {meta_data_importances.sum(axis=1)}\nfragments: {total_fragment_importances.sum(axis=1)}")

# df_feature_importances.to_csv(f'../data/feature_importances_{catalysts[i]}.csv', index=False)








## FROM INITIAL MODELING
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=58)
#
# RF = RandomForestRegressor(n_estimators=1000, max_depth=None, min_samples_split=3, min_samples_leaf=5, max_features=None)
# RF.fit(X_train, y_train)
# y_pred = RF.predict(X_test)
#
# mse = mean_squared_error(y_test, y_pred)
# r_sqr = r2_score(y_test, y_pred)
# print(mse, np.sqrt(mse))
# print(r_sqr)
# train_pred = RF.predict(X_train)
# print(r2_score(y_train, train_pred))
#
# fig, ax = plt.subplots()
# ax.scatter(y_pred, y_test)
# ax.set_xlabel('normalized MALDI intensity')
# ax.set_ylabel('EIC product area')
# ax.set_title('Whole Molecule Noramlzied MALDI vs EIC')
# plt.show()