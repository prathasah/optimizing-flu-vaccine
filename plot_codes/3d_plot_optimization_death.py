import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

df = pd.read_csv("../data/cleaned_optimization_death.csv")
print(df.head())

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(111, projection='3d')

age_groups = df.age_group.unique()
efficacy_groups = df.efficacy.unique()

n=17
my_colors = mpl.cm.get_cmap("Spectral") # Qualitative colormap
color_list = [my_colors(1.*i/n) for i in range(n)]

my_distinct_colors = ["#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe" ,"#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000","#aaffc3","#808000"]

for (num,col) in zip(xrange(n), my_distinct_colors):
	di = df.loc[df['age_group_num'] == num+1]
	xpos = di["age_group_num"]
	ypos = di["efficacy"]
	num_elements = len(xpos)
	zpos = np.zeros(num_elements)
	dx = np.ones(num_elements) *0.5
	dy = np.ones(num_elements)*5
	dz = di["doses"]
	ax1.bar3d(xpos, ypos, zpos, dx,dy,dz, color=col, alpha=0.9)

ax1.set_title('Minimizing mortality', fontsize=22)
ax1.set_xlabel('Age groups', labelpad =10, fontsize=18)
ax1.set_ylabel('Vaccine efficacy (%)', labelpad =10, fontsize=18)
ax1.set_xticks(np.arange(1,19,1))
ax1.set_xticklabels([age_groups[num] for num in np.arange(0,17,1)])
ax1.set_yticks(np.arange(12.5,92.5,10))
ax1.set_yticklabels([efficacy_groups[num] for num in np.arange(0,8,1)])
ax1.set_zlabel('Doses (per 100,000)', fontsize=18)
#plt.show()
ax1.view_init(elev=36., azim=297)
plt.savefig("minimizing_mortality.png",  pad_inches=0.01)
#ax1.view_init(elev=50., azim=320)

