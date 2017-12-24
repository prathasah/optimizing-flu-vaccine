import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

dfi = {}
dfm ={}

dfi  = pd.read_csv("../data/vaccination_results/full_model_dec19_100iter/combined_incidence_11Dec2017.csv")
dfi.relative_coverage = dfi.relative_coverage*100
dfi.vaccine_efficacy = dfi.vaccine_efficacy*100
dfi["infection_per_million"] = dfi.total_infections/1000000
df1 = dfi.loc[(dfi.vaccine_efficacy > 9)]
df1 = df1.loc[(dfi.relative_coverage <100)]
df1.drop_duplicates(subset=['relative_coverage','vaccine_efficacy'], inplace=True)
df1.sort_values(['vaccine_efficacy','relative_coverage'], ascending=[True, True], inplace=True)
df1.to_csv("check_infection.csv")

dfm  =  pd.read_csv("../data/vaccination_results/full_model_dec19_100iter/combined_mortality_11Dec2017.csv")
dfm.relative_coverage = dfm.relative_coverage*100
dfm.vaccine_efficacy = dfm.vaccine_efficacy*100
dfm["mortality_per_thousands"] = dfm.total_mortality/1000
df2 = dfm.loc[(dfm.vaccine_efficacy > 9)]
df2 = df2.loc[(dfm.relative_coverage <100)]
df2.drop_duplicates(subset=['relative_coverage','vaccine_efficacy'], inplace=True)
df2.sort_values(['vaccine_efficacy','relative_coverage'], ascending=[True, True], inplace=True)
df2.to_csv("check_mortality.csv")
#############
#unique elements in x-axis

#print x_elem, y_elem
#xpos = dfi['vaccine_efficacy']
#print ("shape"), xpos
#ypos = dfi['relative_coverage'].as_matrix().reshape(x_elem, y_elem)
#zpos = dfi['infection_per_million'].as_matrix().reshape(x_elem, y_elem)



fig = plt.figure(figsize=(19,7))

ax1 = fig.add_subplot(121, projection='3d')

fig1 = ax1.plot_trisurf(df1.relative_coverage, df1.vaccine_efficacy, df1.infection_per_million, cmap=cm.magma, linewidth=0.2)
ax1.set_xlabel('Population coverage (%)', labelpad =15, fontsize=18)
ax1.set_ylabel('Vaccine efficacy (%)', labelpad =15, fontsize=18)
ax1.set_zlabel('Infections (millions)',labelpad =15, fontsize=18)
ax1.set_xlim(0 ,100)
ax1.set_ylim(10 ,60)
fig.colorbar(fig1, ax=ax1,shrink=0.5, aspect=7, pad =0.09)
ax1.view_init(elev=50., azim=310)


ax2 = fig.add_subplot(122, projection='3d')

fig2 = ax2.plot_trisurf(df2.relative_coverage,df2.vaccine_efficacy, df2.mortality_per_thousands, cmap=cm.viridis, linewidth=0.2)
ax2.set_xlabel('Population coverage (%)', labelpad =15, fontsize=18)
ax2.set_ylabel('Vaccine efficacy (%)', labelpad =15, fontsize=18)
ax2.set_zlabel('Mortality (thousands)',labelpad =15, fontsize=18)
ax2.set_xlim(0 ,100)
ax2.set_ylim(10 ,60)
ax2.set_zlim(0,150)
fig.colorbar(fig2, ax=ax2,shrink=0.5, aspect=7, pad =0.09)
ax2.view_init(elev=50., azim=310)
plt.tight_layout(pad=2, w_pad=0.3, h_pad=0.5)
plt.savefig("../plots/effect_coverage_efficacy_on_burden_100iter.png",  pad_inches=0.0)
plt.savefig("../plots/effect_coverage_efficacy_on_burden_100iter.pdf",  pad_inches=0.0)
