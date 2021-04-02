import pandas as pd
import numpy as np

eval_df = pd.read_csv('../data/predictions.csv')


def rank(m1, m2, m3 , m4):
    array = np.array([m1, m2, m3, m4])
    sorted = np.argsort(array)
    return pd.Series(sorted)
eval_df[['norm_rank1', 'norm_rank2', 'norm_rank3', 'norm_rank4']] = eval_df.apply(lambda x: rank(x['Ru Normalized MALDI'], x['Pd Normalized MALDI'], x['Ir Normalized MALDI'], x['Cu Normalized MALDI']), axis=1)
eval_df[['maldi_rank1', 'maldi_rank2', 'maldi_rank3', 'maldi_rank4']] = eval_df.apply(lambda x: rank(x['Ru MALDI Product Intensity'], x['Pd MALDI Product Intensity'], x['Ir MALDI Product Intensity'], x['Cu MALDI Product Intensity']), axis=1)
eval_df[['true_rank1', 'true_rank2', 'true_rank3', 'true_rank4']] = eval_df.apply(lambda x: rank(x['Ru EIC(+)[M+H] Product Area'], x['Pd EIC(+)[M+H] Product Area'], x['Ir EIC(+)[M+H] Product Area'], x['Cu EIC(+)[M+H] Product Area']), axis=1)
eval_df[['pred_rank1', 'pred_rank2', 'pred_rank3', 'pred_rank4']] = eval_df.apply(lambda x: rank(x['Ru preds'], x['Pd preds'], x['Ir preds'], x['Cu preds']), axis=1)


# preds
eval_df['is_prank1_correct'] = eval_df['true_rank1'] == eval_df['pred_rank1']
eval_df['is_prank2_correct'] = eval_df['true_rank2'] == eval_df['pred_rank2']
eval_df['is_prank3_correct'] = eval_df['true_rank3'] == eval_df['pred_rank3']
eval_df['is_prank4_correct'] = eval_df['true_rank4'] == eval_df['pred_rank4']

eval_df['is_mrank1_correct'] = eval_df['true_rank1'] == eval_df['maldi_rank1']
eval_df['is_mrank2_correct'] = eval_df['true_rank2'] == eval_df['maldi_rank2']
eval_df['is_mrank3_correct'] = eval_df['true_rank3'] == eval_df['maldi_rank3']
eval_df['is_mrank4_correct'] = eval_df['true_rank4'] == eval_df['maldi_rank4']


# pred rankings
print(eval_df['is_prank1_correct'].mean())
print(eval_df['is_prank2_correct'].mean())
print(eval_df['is_prank3_correct'].mean())
print(eval_df['is_prank4_correct'].mean())

print('\n')

# maldi rankings
print(eval_df['is_mrank1_correct'].mean())
print(eval_df['is_mrank2_correct'].mean())
print(eval_df['is_mrank3_correct'].mean())
print(eval_df['is_mrank4_correct'].mean())


eval_df.to_csv('../data/catalyst_rankings.csv', index=False)

