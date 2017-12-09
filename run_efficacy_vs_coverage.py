import efficacy_vs_coverage as efc
import numpy as np
import csv
import random

sys_random = random.SystemRandom()
###########################################3
niter = 5
hosp = np.zeros((niter,17))
inf = np.zeros((niter,17))

##################################################################
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

#for vacEfficacy in np.arange(0.001, 0.61, 0.01):
for vacEfficacy in [0.1]:
	#for relative_coverage in np.arange(0,1, 0.01):
	for relative_coverage in [0.4]:
		for num in xrange(niter):
			totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated, infections, hospitalizations,check =efc.run_efficacy_simulation(relative_coverage, vacEfficacy)
			print ("check!!"), num,check
		
		"""
		elem1 = [relative_coverage,vacEfficacy_2017, sum(list(avg_inf))] + list(avg_inf)
		#print elem1
		writer1.writerow(elem1)
		avg_hosp = np.mean(hosp, axis=0)
		elem2 = [relative_coverage,vacEfficacy_2017, sum(list(avg_hosp))]+list(avg_hosp)
		writer2.writerow(elem2)
		elem3 = [relative_coverage,vacEfficacy_2017, totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated]
		writer3.writerow(elem3)
		"""
