#!/usr/bin/python
import Influenza
import numpy as np
import csv
#########################################
niter=2
##################################################################
#write csv
writer1 = csv.writer(open('Incidence_8Dec2017.csv','wb'))
header = ["relative_coverage","vaccine_efficacy", "total_infections", "age0", "age1", "age5", "age10", "age15", "age20", "age25", "age30", "age35", "age40", "age45", "age50", "age55", "age60", "age65", "age70", "age75"]
writer1.writerow(header)
############################################################
writer2 = csv.writer(open('Hospitalizations_8Dec2017.csv','wb'))
header = ["relative_coverage", "vaccine_efficacy", "total_hospitalizations", "age0", "age1", "age5", "age10", "age15", "age20", "age25", "age30", "age35", "age40", "age45", "age50", "age55", "age60", "age65", "age70", "age75"]
writer2.writerow(header)
############################################
s = Influenza.Simulation()

## typical vaccination proportions: [0.2625729893,0.2016073335,0.1210523011,0.1689583726	0.2458090036]
##for age groups: 6m-4y, 5y - 17y,18y-49y,50y - 64y,65y+
typical_coverage = [0.2625729893,0.2016073335,0.1210523011,0.1689583726, 0.2458090036]

##Vaccine efficacy

## create empty array for hospitlization and infection
hosp = np.zeros((niter,17))
inf = np.zeros((niter,17))
for vacEfficacy_2017 in np.arange(0.001, 0.61, 0.01):
	for relative_coverage in np.arange(0.1,2.1, 0.1):
		coverage_2017 = [relative_coverage*num for num in typical_coverage]
		print ("relative coverage, efficacy"), relative_coverage, vacEfficacy_2017
		for num in xrange(niter):
			## input parameters 1) vac times, 2) proportion vaccinated across age groups, and vaccine efficact
			vacsUsed = s.simulateWithVaccine([0],coverage_2017,vacEfficacy_2017)
			infections, hospitalizations = s.short_output()
			hosp[num] = hospitalizations
			inf[num] = infections	
		avg_inf = np.mean(inf, axis = 0)
		elem1 = [relative_coverage,vacEfficacy_2017, sum(list(avg_inf))] + list(avg_inf)
		#print elem1
		writer1.writerow(elem1)
		avg_hosp = np.mean(hosp, axis=0)
		elem2 = [relative_coverage,vacEfficacy_2017, sum(list(avg_hosp))]+list(avg_hosp)
		writer2.writerow(elem2)

	

