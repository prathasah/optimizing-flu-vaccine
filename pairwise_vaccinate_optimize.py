#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
sys.path.insert(0, './Influenza/Parameters')
import numpy as np
import Simulation_for_pairwise as Simulation
import random
import subprocess

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


def efficacy_optimization_simulation(vacEfficacy, sub_iter):
	vacNumber = 140e6
	vacCoverage = 0.43
	vaccineCoverage = [(vacCoverage*num/np.mean(relative_coverage)) for num in relative_coverage]
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy}, index = sub_iter)
	no_vax_vacsUsed = s.simulateWithVaccine([0], [0]*len(relative_coverage), vacEfficacy)
	no_vax_infections, no_vax_hospitalizations, no_vax_mortality, no_vax_DALY = s.short_output()
	vacsUsed = s.simulateWithVaccine([0], vaccineCoverage, vacEfficacy)
	infections, hospitalizations, mortality,DALY = s.short_output()
	optim_infections = [s.Optimization("totalInfections", vacEfficacy, vacNumber)]
	optim_hospitalizations = [s.Optimization("totalHospitalizations", vacEfficacy, vacNumber)]
	optim_mortality = [s.Optimization("totalDeaths", vacEfficacy, vacNumber)]
	optim_DALY = [s.Optimization("totalDALY", vacEfficacy, vacNumber)]
	return sum(no_vax_infections),  sum(infections), optim_infections,  sum(no_vax_hospitalizations), sum(hospitalizations), optim_hospitalizations, sum(no_vax_mortality),  sum(mortality), optim_mortality, sum(no_vax_DALY),  sum(DALY), optim_DALY
	
		
		
######################################################################################		
if __name__ == "__main__":
    
    
	vacEfficacy_list  = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
	iter_list = [100, 300, 500, 700, 900]
	efficacy_iter_list = []
	for efficacy in vacEfficacy_list:
		for iter1 in iter_list:
			efficacy_iter_list.append((efficacy,iter1))
				


	for index, (vacEfficacy, sub_iter) in enumerate(efficacy_iter_list):
		infections_1, infections_2, optim_infections, hosp_1, hosp_2, optim_hosp, mortality_1, mortality_2, optim_mortality, DALY_1, DALY_2, optim_DALY= efficacy_optimization_simulation(vacEfficacy,1)
		
	
		
