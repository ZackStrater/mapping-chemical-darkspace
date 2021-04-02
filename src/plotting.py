
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')


# baseline RF scores (only using MALDI output)
x = [0, 1, 2, 3]
scores = [0.59, 0.62, 0.71, 0.68]
labels = ['Cu', 'Ir', 'Pd', 'Ru']

fig, ax = plt.subplots()

bars = ax.bar(x, scores)
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
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


fig, ax = plt.subplots()
x = np.array([0, 1, 2, 3, 4])
# 0: just maldi, 1: with metadata, 2: with fragment analysis, 3: with hyperparameter tuning, 4: removing ambiguous structures
Cu_scores = [0.563, 0.71, 0.725, 0.721, 0.723]
Ir_scores = [0.583, 0.683, 0.683, 0.713, 0.739]
Pd_scores = [0.711, 0.791, 0.791, 0.800, 0.818]
Ru_scores = [0.669, 0.763, 0.763, 0.762, 0.78]

ax.plot(x, Cu_scores)
ax.plot(x, Ir_scores)
ax.plot(x, Pd_scores)
ax.plot(x, Ru_scores)
plt.show()


x = np.array([0, 1, 2])

Cu_importances = [0.696, 0.258, 0.046]
Ir_importances = [0.769, 0.180, 0.051]
Pd_importances = [0.804, 0.170, 0.026]
Ru_importances = [0.806, 0.158,  0.035]
features = ['MALDI', 'Metadata', 'Fragments']
ax = plt.subplot(111)
ax.bar(x-0.2, Cu_importances, width=0.2, color='tab:blue', align='center', label='Cu')
ax.bar(x, Ir_importances, width=0.2, color='tab:orange', align='center', label='Ir')
ax.bar(x+0.2, Pd_importances, width=0.2, color='tab:green', align='center', label='Pd')
ax.bar(x+0.4, Ru_importances, width=0.2, color='tab:red', align='center', label='Ru')
ax.set_ylim(0, 1)
ax.set_xticks(x)
ax.set_xticklabels(features, size=15)
ax.set_title('RF Grouped Feature Importances', size=20)
ax.legend()
plt.show()



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
ax.scatter(Ru_x, Ru_y, label='Ru', alpha=0.5, s=20)
ax.scatter(Pd_x, Pd_y, label='Ru', alpha=0.5, s=20)
ax.scatter(Ir_x, Ir_y, label='Ir', alpha=0.5, s=20)
ax.scatter(Cu_x, Cu_y, label='Cu', alpha=0.5, s=20)
ax.set_xlabel('Model Prediction', size=15)
ax.set_ylabel('UPLC-MS EIC-product area', size=15)
ax.set_title('RF Predictions vs UPLC EIC', size=20)
ax.legend()
plt.show()