import numpy as np
import csv
import random
import efficacy_vs_coverage 
import os
###########################################3
niter = 5
hosp = np.zeros((niter,17))
inf = np.zeros((niter,17))
############################################
def delete_cache():

	import os

	# Produces a sorted list of all files in a directory
	dirList = os.listdir("./Influenza")   # Use os.listdir() if want current directory
	dirList.sort()

	# Takes advantage of fact that both py and pyc files will share same name and
	# that pyc files will appear immediately after their py counterparts in dirList

	lastPyName = ""

	for file in dirList:
	    if file[-4:] == ".pyc":
		os.remove("./Influenza/"+str(file))
	
	#print ("check!!"), dirList	
	    #if file[-3:] == ".py":
	    #	lastPyName = file[:-3]
	    #elif file[-4:] == ".pyc":
	    #	if lastPyName == file[:-4]:
		    #os.remove(lastPyName + ".py")
		    #os.remove(lastPyName + ".pyc") # In case you want to delete this too

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
			reload(efficacy_vs_coverage)
			totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated, infections, hospitalizations,check =efficacy_vs_coverage.run_efficacy_simulation(relative_coverage, vacEfficacy)
			print ("check!!"), num,check
			delete_cache()
			
delete_cache()
dirList = os.listdir("./Influenza") 
print dirList		
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
