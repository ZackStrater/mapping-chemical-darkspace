



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