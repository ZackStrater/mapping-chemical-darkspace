

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


xls = pd.ExcelFile('../data/merck_data.xlsx')
df = pd.read_excel(xls, sheet_name=1)

X = df.loc[:, 'Location':'Piperazines N-Aryl']
y_Cu = df.loc[:, 'Cu TWC Product Area': 'Cu Normalized MALDI Pdt/IS']
y_Ir = df.loc[:, 'Ir TWC Product Area': 'Ir Normalized MALDI Pdt/IS']
y_Pd = df.loc[:, 'Pd TWC Product Area': 'Pd Normalized MALDI Pdt/IS']
y_Ru = df.loc[:, 'Ru TWC Product Area': 'Ru Normalized MALDI Pdt/IS']

blank_cat_df = pd.DataFrame(np.zeros((384, 4)), columns=['Cu_cat', 'Ir_cat', 'Pd_cat', 'Ru_cat'])

y_Cu = pd.concat([y_Cu, blank_cat_df], axis=1)
y_Cu['Cu_cat'] = 1
y_Ir = pd.concat([y_Ir, blank_cat_df], axis=1)
y_Ir['Ir_cat'] = 1
y_Pd = pd.concat([y_Pd, blank_cat_df], axis=1)
y_Pd['Pd_cat'] = 1
y_Ru = pd.concat([y_Ru, blank_cat_df], axis=1)
y_Ru['Ru_cat'] = 1

df_Cu = pd.concat([X, y_Cu], axis=1)
df_Cu.rename(columns=dict(zip(df_Cu.columns[-12:-4], df_Cu.columns[-12:-4].str[3:])), inplace=True)
df_Ir = pd.concat([X, y_Ir], axis=1)
df_Ir.rename(columns=dict(zip(df_Ir.columns[-12:-4], df_Ir.columns[-12:-4].str[3:])), inplace=True)
df_Pd = pd.concat([X, y_Pd], axis=1)
df_Pd.rename(columns=dict(zip(df_Pd.columns[-12:-4], df_Pd.columns[-12:-4].str[3:])), inplace=True)
df_Ru = pd.concat([X, y_Ru], axis=1)
df_Ru.rename(columns=dict(zip(df_Ru.columns[-12:-4], df_Ru.columns[-12:-4].str[3:])), inplace=True)

df_out = pd.concat([df_Cu, df_Ir, df_Pd, df_Ru], axis=0)
print(df_out)

df_out.to_csv('../data/merck_data_combined', index=False)