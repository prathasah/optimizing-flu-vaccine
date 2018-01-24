import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
##################################################
def get_axis_limits(ax, scale=.7):
    return ax.get_xlim()[0]*scale, ax.get_ylim()[1]*scale, ax.get_zlim()[1]*scale
#################################################################3

dfi  = pd.read_csv("../data/vaccination_results/combined_incidence_5Jan2018_500iter.csv")
dfi.relative_coverage = dfi.relative_coverage*100
dfi.vaccine_efficacy = dfi.vaccine_efficacy*100
dfi['infection_per_million'] = dfi.total_infections/1000000
df1 = dfi.loc[(dfi.vaccine_efficacy > 9)]
df1 = df1.loc[(dfi.relative_coverage <100)]
df1.drop_duplicates(subset=['relative_coverage','vaccine_efficacy'], inplace=True)
df1.sort_values(['vaccine_efficacy','relative_coverage'], ascending=[True, True], inplace=True)
print max(dfi['infection_per_million'])
#df1.to_csv("check_infection.csv")


dfh  =  pd.read_csv("../data/vaccination_results/combined_hospitalization_5Jan2018_500iter.csv")
dfh.relative_coverage = dfh.relative_coverage*100
dfh.vaccine_efficacy = dfh.vaccine_efficacy*100
dfh["hosp_per_million"] = dfh.total_hospitalizations/1000000
df3 = dfh.loc[(dfh.vaccine_efficacy > 9)]
df3 = df3.loc[(dfh.relative_coverage <100)]
df3.drop_duplicates(subset=['relative_coverage','vaccine_efficacy'], inplace=True)
df3.sort_values(['vaccine_efficacy','relative_coverage'], ascending=[True, True], inplace=True)
#df3.to_csv("check_hospitalization.csv")

dfm  =  pd.read_csv("../data/vaccination_results/combined_mortality_5Jan2018_500iter.csv")
dfm.relative_coverage = dfm.relative_coverage*100
dfm.vaccine_efficacy = dfm.vaccine_efficacy*100
dfm["mortality_per_thousands"] = dfm.total_mortality/1000
df2 = dfm.loc[(dfm.vaccine_efficacy > 9)]
df2 = df2.loc[(dfm.relative_coverage <100)]
df2.drop_duplicates(subset=['relative_coverage','vaccine_efficacy'], inplace=True)
df2.sort_values(['vaccine_efficacy','relative_coverage'], ascending=[True, True], inplace=True)
#df2.to_csv("check_mortality.csv")


dfd  =  pd.read_csv("../data/vaccination_results/combined_DALY_5Jan2018_500iter.csv")
dfd.relative_coverage = dfd.relative_coverage*100
dfd.vaccine_efficacy = dfd.vaccine_efficacy*100
dfd["DALY_per_million"] = dfd.total_DALY/1000000
df4 = dfd.loc[(dfd.vaccine_efficacy > 9)]
df4 = df4.loc[(dfd.relative_coverage <100)]
df4.drop_duplicates(subset=['relative_coverage','vaccine_efficacy'], inplace=True)
df4.sort_values(['vaccine_efficacy','relative_coverage'], ascending=[True, True], inplace=True)
df4.to_csv("check_DALY.csv")

#############
#unique elements in x-axis

#print x_elem, y_elem
#xpos = dfi['vaccine_efficacy']
#print ("shape"), xpos
#ypos = dfi['relative_coverage'].as_matrix().reshape(x_elem, y_elem)
#zpos = dfi['infection_per_million'].as_matrix().reshape(x_elem, y_elem)



fig = plt.figure(figsize=(20,14))

ax1 = fig.add_subplot(221, projection='3d')

print min(df4.DALY_per_million)
print max(df4.DALY_per_million)
fig1 = ax1.plot_trisurf(df1.relative_coverage, df1.vaccine_efficacy, df1.infection_per_million, cmap=cm.magma, linewidth=0.2)
ax1.set_xlabel('Population coverage (%)', labelpad =15, fontsize=18)
ax1.set_ylabel('Vaccine efficacy (%)', labelpad =15, fontsize=18)
ax1.set_zlabel('Incidence (millions)',labelpad =15, fontsize=18)
ax1.set_xlim(0 ,100)
ax1.set_ylim(10 ,80)
fig.colorbar(fig1, ax=ax1,shrink=0.5, aspect=7, pad =0.09)
ax1.view_init(elev=50., azim=310)


ax2 = fig.add_subplot(222, projection='3d')

fig2 = ax2.plot_trisurf(df3.relative_coverage, df3.vaccine_efficacy, df3.hosp_per_million, cmap=cm.RdPu_r, linewidth=0.2)
ax2.set_xlabel('Population coverage (%)', labelpad =15, fontsize=18)
ax2.set_ylabel('Vaccine efficacy (%)', labelpad =15, fontsize=18)
ax2.set_zlabel('Hospitalizations (millions)',labelpad =15, fontsize=18)
ax2.set_xlim(0 ,100)
ax2.set_ylim(10 ,80)
fig.colorbar(fig2, ax=ax2,shrink=0.5, aspect=7, pad =0.09)
ax2.view_init(elev=50., azim=310)


ax3 = fig.add_subplot(223, projection='3d')

fig3 = ax3.plot_trisurf(df2.relative_coverage,df2.vaccine_efficacy, df2.mortality_per_thousands, cmap=cm.viridis, linewidth=0.2)
ax3.set_xlabel('Population coverage (%)', labelpad =15, fontsize=18)
ax3.set_ylabel('Vaccine efficacy (%)', labelpad =15, fontsize=18)
ax3.set_zlabel('Mortality (thousands)',labelpad =15, fontsize=18)
ax3.set_xlim(0 ,100)
ax3.set_ylim(10 ,80)
ax3.set_zlim(0,150)
fig.colorbar(fig3, ax=ax3,shrink=0.5, aspect=7, pad =0.09)
ax3.view_init(elev=50., azim=310)
plt.tight_layout(pad=4, w_pad=0.3, h_pad=0.5)

ax4 = fig.add_subplot(224, projection='3d')
fig4 = ax4.plot_trisurf(df4.relative_coverage,df4.vaccine_efficacy, df4.DALY_per_million, cmap=cm.BuPu_r, linewidth=0.2)
ax4.set_xlabel('Population coverage (%)', labelpad =15, fontsize=18)
ax4.set_ylabel('Vaccine efficacy (%)', labelpad =15, fontsize=18)
ax4.set_zlabel('DALYs (millions)',labelpad =15, fontsize=18)
ax4.set_xlim(0 ,100)
ax4.set_ylim(10, 80)
ax4.set_zlim(0,7.5)
fig.colorbar(fig4, ax=ax4,shrink=0.5, aspect=7, pad =0.09)
ax4.view_init(elev=50., azim=310)


plt.tight_layout(pad=3, w_pad=0.3, h_pad=2)
print get_axis_limits(ax1)
ax1.text(-80, 50, 100, "A", None, size =20)
ax2.text(-80, 50, 1.3, "B", None, size =20)
ax3.text(-80, 50, 80, "C", None, size =20)
ax4.text(-80, 50, 40, "D", None, size =20)
plt.savefig("../plots/effect_typical_coverage_efficacy_on_outcome.png",  pad_inches=0.0)
plt.savefig("../plots/effect_typical_coverage_efficacy_on_outcome.pdf",  pad_inches=0.0, dpi= 100)
#plt.savefig("../plots/effect_typical_coverage_efficacy_on_outcome.eps",  pad_inches=0.0, dpi=1000)
