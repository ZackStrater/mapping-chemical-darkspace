

import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

xls = pd.ExcelFile('../data/merck_data.xlsx')
df = pd.read_excel(xls, sheet_name=1)
print(df.head())

