import numpy as np
import csv
import random
import efficacy_vs_coverage 
import os
###########################################3
niter = 10

############################################
#write csv
writer1 = csv.writer(open('Incidence_8Dec2017.csv','wb'))
header = ["relative_coverage","vaccine_efficacy", "total_infections", "age0", "age1", "age5", "age10", "age15", "age20", "age25", "age30", "age35", "age40", "age45", "age50", "age55", "age60", "age65", "age70", "age75"]
writer1.writerow(header)
############################################################
writer2 = csv.writer(open('Hospitalizations_8Dec2017.csv','wb'))
header = ["relative_coverage", "vaccine_efficacy", "total_hospitalizations", "age0", "age1", "age5", "age10", "age15", "age20", "age25", "age30", "age35", "age40", "age45", "age50", "age55", "age60", "age65", "age70", "age75"]
writer2.writerow(header)

##################################################
writer3 = csv.writer(open('vaccinated_comparisons_8Dec2017.csv','wb'))
header = ["relative_coverage", "vaccine_efficacy", "total_population","total_unvaccinated", "total_vaccinated", "unvaccinated_total_infections", "vaccinated_total_infections"]
writer3.writerow(header)
############################################

for vacEfficacy in np.arange(0.001, 0.61, 0.01):
	for relative_coverage in np.arange(0,1, 0.01):
		hosp = np.zeros((niter,17))
		inf = np.zeros((niter,17))
		totpop = []
		tot_vax = []
		tot_unvax = []
		vax_inf = []
		unvax_inf = []	
		print vacEfficacy, relative_coverage
		for num in xrange(niter):
			tot_population, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated, infections, hospitalizations,check =efficacy_vs_coverage.run_efficacy_simulation(relative_coverage, vacEfficacy)
			hosp[num] = hospitalizations
			inf[num] = infections
			totpop.append(tot_population)
			tot_vax.append(tot_vaccinated)
			tot_unvax.append(tot_unvaccinated)
			vax_inf.append(vaccinated)
			unvax_inf.append(unvaccinated)
			
			
		avg_inf = np.mean(inf, axis = 0)
		avg_hosp = np.mean(hosp, axis=0)
		elem1 = [relative_coverage,vacEfficacy, sum(list(avg_inf))] + list(avg_inf)
		writer1.writerow(elem1)
		avg_hosp = np.mean(hosp, axis=0)
		elem2 = [relative_coverage,vacEfficacy, sum(list(avg_hosp))]+list(avg_hosp)
		writer2.writerow(elem2)
		elem3 = [relative_coverage,vacEfficacy, np.mean(totpop), np.mean(tot_unvax), np.mean(tot_vax), np.mean(unvax_inf), np.mean(vax_inf)]
		writer3.writerow(elem3)

