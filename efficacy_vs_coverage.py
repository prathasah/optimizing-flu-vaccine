#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import numpy as np
import Simulation 


def run_efficacy_simulation(coverage, vacEfficacy):
	reload(Simulation)
	s = Simulation.run_Simulation(paramValues = {"vacEfficacy":vacEfficacy, "vacCoverage": coverage})
	coverage = s.vaccine_coverage()
	vacsUsed = s.simulateWithVaccine([0], coverage, vacEfficacy)
	infections, hospitalizations, mortality = s.short_output()
	check = s.debug_info()
	totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated = s.vaccinated_output()
	return totpop, tot_unvaccinated, tot_vaccinated, unvaccinated, vaccinated, infections, hospitalizations, mortality, check
	
		
		
	

