
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')


# baseline RF scores (only using MALDI output)
x = [0, 1, 2, 3]
scores = [0.56, 0.58, 0.71, 0.67]
labels = ['Cu', 'Ir', 'Pd', 'Ru']

fig, ax = plt.subplots()

bars = ax.bar(x, scores)
colors = ['tab:green', 'tab:blue', 'tab:red', 'tab:orange']
for i in range(4):
    bars[i].set_color(colors[i])
for i, v in enumerate(scores):
    ax.text(x[i]-0.2, v + 0.01, str(v), size=15)

ax.set_xticks(x)
ax.set_xticklabels(labels, size=15)
ax.set_ylim(0, 1)
ax.set_ylabel('R^2 score', size=15)
ax.set_title('Random Forest with MALDI Output', size=20)
plt.show()


# model improvement
fig, ax = plt.subplots(figsize=(10, 8))
x = np.array([0, 1, 2, 3, 4])
# 0: just maldi, 1: with metadata, 2: with fragment analysis, 3: with hyperparameter tuning, 4: removing ambiguous structures
Cu_scores = [0.563, 0.71, 0.725, 0.721, 0.723]
Ir_scores = [0.583, 0.683, 0.683, 0.713, 0.739]
Pd_scores = [0.711, 0.791, 0.791, 0.800, 0.818]
Ru_scores = [0.669, 0.763, 0.763, 0.762, 0.78]
ax.plot(x, Cu_scores, color='tab:green', label='Cu R^2')
ax.plot(x, Ir_scores, color='tab:blue', label='Ir R^2')
ax.plot(x, Pd_scores, color='tab:red', label='Pd R^2')
ax.plot(x, Ru_scores, color='tab:orange', label='Ru R^2')
ax.scatter(x, Cu_scores, color='tab:green')
ax.scatter(x, Ir_scores, color='tab:blue')
ax.scatter(x, Pd_scores, color='tab:red')
ax.scatter(x, Ru_scores, color='tab:orange')
ax.set_ylim(0.55, 0.95)
ax.axhline(y=0.85, color='red', ls='--')
ax.set_xlim(-0.5, 4.5)
ax.set_xticks(x)
ax.set_xticklabels(['MALDI\nonly', '+ Molecule\nMetadata', '+ Fragment\nAnalysis', '+ Hyperparameter\nTuning', '+  Removing\nAmbiguous\nStructures'], size=13)
ax.text(2.02, 0.86, 'approximate accuracy limit\ndue to irreducible error', size=15)
ax.set_ylabel('RF R^2 Score', size=20)
ax.set_title('RF Model Improvement', size=25)
ax.text(-0.3, 0.56, '0.56', size=12)
ax.text(-0.3, 0.58, '0.58', size=12)
ax.text(-0.3, 0.708, '0.71', size=12)
ax.text(-0.3, 0.666, '0.67', size=12)
ax.text(4.1, 0.72, '0.72', size=12)
ax.text(4.1, 0.739, '0.74', size=12)
ax.text(4.1, 0.818, '0.82', size=12)
ax.text(4.1, 0.780, '0.78', size=12)

ax.legend()
plt.show()




#plotting feature importances
x = np.array([0, 1, 2])

Cu_importances = [0.696, 0.258, 0.046]
Ir_importances = [0.769, 0.180, 0.051]
Pd_importances = [0.804, 0.170, 0.026]
Ru_importances = [0.806, 0.158,  0.035]
features = ['MALDI\nOutput', 'Molecule\nMetadata', 'Fragments\nAnalysis']
ax = plt.subplot(111)
ax.bar(x-0.2, Cu_importances, width=0.2, color='tab:green', align='center', label='Cu')
ax.bar(x, Ir_importances, width=0.2, color='tab:blue', align='center', label='Ir')
ax.bar(x+0.2, Pd_importances, width=0.2, color='tab:red', align='center', label='Pd')
ax.bar(x+0.4, Ru_importances, width=0.2, color='tab:orange', align='center', label='Ru')
ax.set_ylim(0, 1)
ax.set_xticks(x)
ax.set_xticklabels(features, size=15)
ax.set_title('RF Grouped Feature Importances', size=20)
ax.set_ylabel('Gini Importances', size=15)
ax.legend()
plt.show()


# predictions
df = pd.read_csv('../data/predictions.csv')
df_Cu = df[df['Cu_cat'] == 1].reset_index(drop=True)
df_Ir = df[df['Ir_cat'] == 1].reset_index(drop=True)
df_Pd = df[df['Pd_cat'] == 1].reset_index(drop=True)
df_Ru = df[df['Ru_cat'] == 1].reset_index(drop=True)
# Cu
Cu_x = df_Cu['preds']
Cu_y = df_Cu['EIC(+)[M+H] Product Area']
Ir_x = df_Ir['preds']
Ir_y = df_Ir['EIC(+)[M+H] Product Area']
Pd_x = df_Pd['preds']
Pd_y = df_Pd['EIC(+)[M+H] Product Area']
Ru_x = df_Ru['preds']
Ru_y = df_Ru['EIC(+)[M+H] Product Area']
pred_EIC, ax = plt.subplots()
ax.scatter(Ru_x, Ru_y, label='Ru', alpha=0.5, s=20, color='tab:orange')
ax.scatter(Pd_x, Pd_y, label='Ru', alpha=0.5, s=20, color='tab:red')
ax.scatter(Ir_x, Ir_y, label='Ir', alpha=0.5, s=20, color='tab:blue')
ax.scatter(Cu_x, Cu_y, label='Cu', alpha=0.5, s=20, color='tab:green')
ax.set_xlabel('Model Prediction', size=15)
ax.set_ylabel('UPLC-MS EIC-Product Area', size=15)
ax.set_ylim(-0.3e6, 7e6)
ax.set_xlim(-0.2e6, 7e6)
ax.plot([0,7e6],[0,7e6], transform=ax.transAxes, color='black', linestyle='--', alpha=0.5)
ax.set_title('RF Predictions vs UPLC EIC', size=20)
ax.legend()
plt.show()