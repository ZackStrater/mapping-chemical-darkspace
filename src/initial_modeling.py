

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


df = pd.read_csv('../data/merck_data_combined.csv')

y = df.pop('EIC(+)[M+H] Product Area')
X = df.copy()

X.drop(['TWC Product Area', 'TWC Product Area%', 'MS TIC(+) Product Area%', 'MS TIC x TWC',
        'Location', 'Sample Name', 'bromide or amine', 'CHEMISTRY', 'Canonical_Smiles'], axis=1, inplace=True)
X['3/4-pyridyl-type 6 MR heterocycles (Br only)'].fillna(0, inplace=True)

X = X[0:384]
y = y[0:384]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=58)

RF = RandomForestRegressor(n_estimators=1000, max_depth=None, min_samples_split=3, min_samples_leaf=5, max_features=None)
RF.fit(X_train, y_train)
y_pred = RF.predict(X_test)

mse = mean_squared_error(y_test, y_pred)
r_sqr = r2_score(y_test, y_pred)
print(mse, np.sqrt(mse))
print(r_sqr)
train_pred = RF.predict(X_train)
print(r2_score(y_train, train_pred))

fig, ax = plt.subplots()
ax.scatter(y_pred, y_test)
ax.set_xlabel('normalized MALDI intensity')
ax.set_ylabel('EIC product area')
ax.set_title('Whole Molecule Noramlzied MALDI vs EIC')
plt.show()