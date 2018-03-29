#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation

# source :  https://www.cdc.gov/flu/fluvaxview/index.htm	
#------------------------------------------------------------
#Year	|  6m-4y    |5y - 17y  | 18y-49y |50y - 64y|	65y+|
#-----------------------------------------------------------
#2012-13|	69.8|	53.1   | 31.1	 |45.1	   |	66.2
#2013-14|	70.4|	55.3   |32.3	 |45.3     |	65
#2014-15|	70.4|	55.8   |33.5	 |47	   |	66.7
#2015-16|	70  |	55.9   |32.7	 |43.6	   | 	63.4
#2016-17|	70  | 	55.6   |33.6	 |45.4	   |	65.3
#----------------------------------------------------------					
					
relative_coverage = [ 0.7012,  0.7012,  0.5514,  0.5514 , 0.4614,  0.3264,  0.3264,  0.3264,  0.3264,
  0.3264,  0.4528,  0.4528,  0.4528,  0.4528,  0.6532,  0.6532,  0.6532]


def run_efficacy_simulation(coverage, vacEfficacy, sub_iter):
	reload(Simulation)
	vaccineCoverage = [(coverage*num/np.mean(relative_coverage)) for num in relative_coverage]
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy}, index = sub_iter)
	vacsUsed = s.simulateWithVaccine([0], vaccineCoverage, vacEfficacy)
	infections, hospitalizations, mortality,DALY = s.short_output()
	#check = s.debug_info()
	#totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated = s.vaccinated_output()
	return infections, hospitalizations, mortality, DALY
	
		
		

######################################################################33
if __name__ == "__main__":
	
	"""
	no_vax_incidence = []
	for sub_index in xrange(1000):
		infections, hospitalizations, mortality, DALY = run_efficacy_simulation(0, 0.6,sub_index)
		#print num
		no_vax_incidence.append(sum(infections)/1e6)
		
	print ("no vaccination incidence"), np.mean(no_vax_incidence), np.std(no_vax_incidence)
	"""
	
	vax_60_incidence = []
	for sub_index in xrange(1):
		#print sub_index
		infections, hospitalizations, mortality, DALY = run_efficacy_simulation(0.43, 1,sub_index)
		vax_60_incidence.append(sum(infections)/1e6)
		
	print ("typical vaccination incidence"), np.mean(vax_60_incidence), np.std(vax_60_incidence)
	