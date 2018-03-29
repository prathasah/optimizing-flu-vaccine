import csv
import numpy


######################################################################################		
if __name__ == "__main__":
    
    
    header = ["iter", "recovery_rate", "latency_rate", "susceptibility_0", "susceptibility_5", "susceptibility_25", "susceptibility_50", "susceptibility_65", "relative_vaccineEfficacyVsInfection_0", "relative_vaccineEfficacyVsInfection_0.5","relative_vaccineEfficacyVsInfection_16","relative_vaccineEfficacyVsInfection_65","vaccineEfficacyVsDeath_0", "vaccineEfficacyVsDeath_20","vaccineEfficacyVsDeath_65", "vaccineEfficacyVsHospitalization_0", "vaccineEfficacyVsHospitalization_5", "vaccineEfficacyVsHospitalization_16", "vaccineEfficacyVsHospitalization_65","caseMortality_0","caseMortality_5","caseMortality_18","caseMortality_50","caseMortality_65",  "caseHospitalization_0","caseHospitalization_5","caseHospitalization_18","caseHospitalization_50","caseHospitalization_65", "R0"]
    writer = csv.writer(open('sampled_parameter_set.csv','wb'))
    writer.writerow(header)
    
    for num in xrange(1000):
        rec_rate = 1. / numpy.random.triangular(0.9, 2.5, 4.6)
        lat_rate = 1. / numpy.random.triangular(1.1, 2.0, 2.5)
        susc_0 = numpy.random.uniform(0.9, 1.)
        susc_5 =      numpy.random.uniform(0.7, 1.)
        susc_25 =     numpy.random.uniform(0.6, 1.)
        susc_50 =     numpy.random.uniform(0.5, 1.)
        susc_65 =      numpy.random.uniform(0.3, 1)
        vac_eff_inf_0 = 0
        vac_eff_inf_0_5 = numpy.random.triangular(0.4, 0.7, 1.)
        vac_eff_inf_16 = numpy.random.triangular(0.5, 0.7, 0.9)
        vac_eff_inf_65 = numpy.random.triangular(0.4, 0.5, 0.6)
        vac_eff_death_0 = numpy.random.uniform(0.4, 0.75)
        vac_eff_death_20 = numpy.random.uniform(0.4, 0.7)
        vac_eff_death_65 = numpy.random.uniform(0.3, 0.7)
        vac_eff_hosp_0 = numpy.random.triangular(0.44, 0.6, 0.72)
        vac_eff_hosp_5 = numpy.random.triangular(0.35, 0.53,0.66)
        vac_eff_hosp_16 = numpy.random.triangular(0.24, 0.43, 0.57)
        vac_eff_hosp_65 = numpy.random.triangular(0.37, 0.50, 0.60)
        case_mort_0 = numpy.random.triangular(0.00002, 0.00004, 0.00026)
        case_mort_5 = numpy.random.triangular(0.00001, 0.00002, 0.00010)
        case_mort_18 = numpy.random.triangular(0.00009, 0.00010, 0.00159)
        case_mort_50 = numpy.random.triangular(0.00010, 0.00134, 0.00159)
        case_mort_65 = numpy.random.triangular(0.00010, 0.00090, 0.01170)
        case_hosp_0 = numpy.random.triangular(0.0033, 0.0141, 0.0245)
        case_hosp_5 = numpy.random.triangular(0.0006, 0.0011, 0.0061)
        case_hosp_18 = numpy.random.triangular(0.0015, 0.0042, 0.0300)
        case_hosp_50 = numpy.random.triangular(0.0015, 0.0193, 0.0300)
        case_hosp_65 = numpy.random.triangular(0.0016, 0.0184, 0.0421)
        R0 = 1.2375
        
        elements = [num, rec_rate, lat_rate, susc_0, susc_5, susc_25, susc_50, susc_65, vac_eff_inf_0, vac_eff_inf_0_5, vac_eff_inf_16, vac_eff_inf_65, vac_eff_death_0, vac_eff_death_20, vac_eff_death_65, vac_eff_hosp_0, vac_eff_hosp_5, vac_eff_hosp_16, vac_eff_hosp_65, case_mort_0, case_mort_5, case_mort_18, case_mort_50, case_hosp_65, case_hosp_0, case_hosp_5, case_hosp_18, case_hosp_50, case_hosp_65, R0]
        writer.writerow(elements)
        
        


        
        
        
