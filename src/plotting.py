
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')

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

