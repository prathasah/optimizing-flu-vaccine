#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation 
## typical vaccination proportions: [0.7012, 0.5514, 0.3264, 0.4528, 0.6532]
##for age groups: 6m-4y, 5y - 17y,18y-49y,50y - 64y,65y+
typical_coverage = [0.7012, 0.5514, 0.3264, 0.4528, 0.6532]


def run_efficacy_simulation(relative_coverage, vacEfficacy):
	reload(Simulation)
	absolute_coverage = [(relative_coverage*num/np.mean(typical_coverage)) for num in typical_coverage]
	s = Simulation.run_Simulation()
	vacsUsed = s.simulateWithVaccine([0],absolute_coverage, vacEfficacy )
	infections, hospitalizations = s.short_output()
	check = s.debug_info()
	totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated = s.vaccinated_output()
	return totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated, infections, hospitalizations, check
		
		
		
	

