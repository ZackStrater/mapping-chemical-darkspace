
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

xls = pd.ExcelFile('../data/merck_data.xlsx')
df = pd.read_excel(xls, sheet_name=1)
df_compare = df[['Ru Normalized MALDI Pdt/IS', 'Pd Normalized MALDI Pdt/IS', 'Ir Normalized MALDI Pdt/IS',
                 'Cu Normalized MALDI Pdt/IS', 'Ru EIC(+)[M+H] Product Area', 'Pd EIC(+)[M+H] Product Area',
                 'Ir EIC(+)[M+H] Product Area', 'Cu EIC(+)[M+H] Product Area']]

def rank(m1, m2, m3 , m4):
    array = np.array([m1, m2, m3, m4])
    sorted = np.argsort(array)
    return pd.Series(sorted)
df_compare[['a', 'b', 'c', 'd']] = df_compare.apply(lambda x: rank(x['Ru Normalized MALDI Pdt/IS'], x['Pd Normalized MALDI Pdt/IS'], x['Ir Normalized MALDI Pdt/IS'], x['Cu Normalized MALDI Pdt/IS']), axis=1)
df_compare[['aa', 'bb', 'cc', 'dd']] = df_compare.apply(lambda x: rank(x['Ru EIC(+)[M+H] Product Area'], x['Pd EIC(+)[M+H] Product Area'], x['Ir EIC(+)[M+H] Product Area'], x['Cu EIC(+)[M+H] Product Area']), axis=1)

df_compare.to_csv('../data/compare.csv', index=False)

# EDA visualization

Ru_MALDI = df['Ru MALDI Product Intensity']
Ru_norm_MALDI = df['Ru Normalized MALDI Pdt/IS']
Ru_EIC = df['Ru EIC(+)[M+H] Product Area']

Pd_MALDI = df['Pd MALDI Product Intensity']
Pd_norm_MALDI = df['Pd Normalized MALDI Pdt/IS']
Pd_EIC = df['Pd EIC(+)[M+H] Product Area']

Ir_MALDI = df['Ir MALDI Product Intensity']
Ir_norm_MALDI = df['Ir Normalized MALDI Pdt/IS']
Ir_EIC = df['Ir EIC(+)[M+H] Product Area']

Cu_MALDI = df['Cu MALDI Product Intensity']
Cu_norm_MALDI = df['Cu Normalized MALDI Pdt/IS']
Cu_EIC = df['Cu EIC(+)[M+H] Product Area']

norm_MALDI_EIC, ax = plt.subplots()
ax.scatter(Ru_norm_MALDI, Ru_EIC, label='Ru', alpha=0.5, s=20)
ax.scatter(Pd_norm_MALDI, Pd_EIC, label='Pd', alpha=0.5, s=20)
ax.scatter(Ir_norm_MALDI, Ir_EIC, label='Ir', alpha=0.5, s=20)
ax.scatter(Cu_norm_MALDI, Cu_EIC, label='Cu', alpha=0.5, s=20)
ax.set_xlabel('normalized MALDI intensity')
ax.set_ylabel('EIC product area')
ax.set_title('normalized MALDI vs EIC')
ax.legend()


MALDI_EIC, ax2 = plt.subplots()
ax2.scatter(Ru_MALDI, Ru_EIC, label='Ru', alpha=0.5, s=20)
ax2.scatter(Pd_MALDI, Pd_EIC, label='Ru', alpha=0.5, s=20)
ax2.scatter(Ir_MALDI, Ir_EIC, label='Ir', alpha=0.5, s=20)
ax2.scatter(Cu_MALDI, Cu_EIC, label='Cu', alpha=0.5, s=20)
ax2.set_xlabel('raw MALDI intensity')
ax2.set_ylabel('EIC product area')
ax2.set_title('raw MALDI vs EIC')
ax2.legend()
plt.show()


# combining data from each catalyst

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

#df_out.to_csv('../data/merck_data_combined.csv', index=False)


