# -*- coding: iso-8859-1 -*-
#
# Epidemiological parameter values
#

from PiecewiseAgeParameter import PiecewiseAgeRate
import pandas as pd
df  = pd.read_csv("/Users/prathasah/Dropbox (Bansal Lab)/Git-files/src/sampled_parameter_set.csv")


def recoveryRatePW(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return PiecewiseAgeRate(
    [df.at[index, "recovery_rate"]],
    [0])


def latencyRatePW(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return PiecewiseAgeRate(
    [df.at[index, "latency_rate"]],
    [0])

def susceptibilityPW(index):
    return PiecewiseAgeRate(
    [df.at[index, "susceptibility_0"],
     df.at[index, "susceptibility_5"],
     df.at[index, "susceptibility_25"],
     df.at[index, "susceptibility_50"],
     df.at[index, "susceptibility_65"]],
    [0, 5, 25, 50, 65])

def  relative_vaccineEfficacyVsInfectionPW(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return PiecewiseAgeRate(
    [df.at[index, "relative_vaccineEfficacyVsInfection_0"],
     df.at[index, "relative_vaccineEfficacyVsInfection_0.5"],
     df.at[index, "relative_vaccineEfficacyVsInfection_16"],
     df.at[index, "relative_vaccineEfficacyVsInfection_65"]],
    [0, 0.5,  16, 65])

def vaccineEfficacyVsDeathPW(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsDeath_0"],
     df.at[index, "vaccineEfficacyVsDeath_20"],
     df.at[index, "vaccineEfficacyVsDeath_65"]],
    [0,20,65])

def vaccineEfficacyVsHospitalizationPW(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return PiecewiseAgeRate(
    [df.at[index, "vaccineEfficacyVsHospitalization_0"],
     df.at[index, "vaccineEfficacyVsHospitalization_5"],
     df.at[index, "vaccineEfficacyVsHospitalization_16"],
     df.at[index, "vaccineEfficacyVsHospitalization_65"]],
    [0, 5, 16, 65])


def caseMortalityPW(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return PiecewiseAgeRate(
    [df.at[index, "caseMortality_0"],
     df.at[index, "caseMortality_5"],
     df.at[index, "caseMortality_18"],
     df.at[index, "caseMortality_50"],
    df.at[index, "caseMortality_65"]],
    [0, 5, 18, 50, 65])


def caseHospitalizationPW(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return PiecewiseAgeRate(
    [df.at[index, "caseHospitalization_0"],
     df.at[index, "caseHospitalization_5"],
     df.at[index, "caseHospitalization_18"],
     df.at[index, "caseHospitalization_50"],
    df.at[index, "caseHospitalization_65"],],
    [0, 5, 18, 50, 65])


def R0(index):
    #df  = pd.read_csv("sampled_parameter_set.csv")
    return df.at[index, "R0"]


