#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:33:38 2024

@author: astronish16
"""
import sys
import json
from icecream import ic
import numpy as np
import scipy.stats as st
import astropy.units as u
from datetime import datetime
import matplotlib.pyplot as plt


import functions as f


# Read the configuration file for DBM/ P-DBM application
file = open('Boundary_condition.json', 'r')
config = json.load(file)

# DBM run type
PDBM = config['P-DBM']
N = config['Number_of_Runs']


# CME Initial Condition
Time_UTC = config['Event_Date']
ic(Time_UTC)
R0 = config['Initial_Position'] * u.R_sun
V0 = config['Initial_Speed'] * u.km/u.s
dR0 = config['Error_Initial_Position'] * u.R_sun
dV0 = config['Error_Initial_Speed']*u.km/u.s
dT0 = config['Possible_Event_delay'] * u.s

# Solar Wind Condition
Auto_DBM = config['Auto_DBM_Parameters']
Wind_type = config['Wind_Type']

# CME cone property
TwoD = config['2D_CONE_DBM']
Omega = config['Half_Width_CME']
CME_Propagation = config['CME_Propagation_Direction']
Phi_CME = f.Phi_Correction(CME_Propagation, Time_UTC)

# Target Information
Target = config['Target']
Phi_target, R_target = f.position(Target, Time_UTC)

# Conversion for Calculation
Event_Time = datetime.strptime(Time_UTC, "%Y-%m-%d %H:%M")
ic(Event_Time)
R0 = (R0.to(u.km)).value
dR0 = (dR0.to(u.km)).value
V0 = V0.value
dV0 = dV0.value
R1 = ((R_target*u.au).to(u.km)).value
T0 = Event_Time.timestamp()
dT0 = dT0.value


# =============================================================================
# Solar wind Properties input
# =============================================================================
'''
First we need to verify the correct solar wind type is passed or not in the boundary condition
If Auto_DBM is True then we need to program should get value of solar wind speed PDF obtained in (Mugatwala et. al, 2024) paper.
In case of P-DBM it will pass the array of solar wind speed other wise mean value of PDF described.
'''

# Verification of correct solar wind type.
if (Wind_type != "Slow") and (Wind_type != "Fast"):
    ic()
    print(f'Wind Type is: {Wind_type}')
    print("Invalid Solar Wind Type")
    print("Please Select 'Slow' or 'Fast' in the Boundary Condition File.")
    sys.exit()


if (PDBM == True) and (Auto_DBM == True):
    ic(PDBM)
    ic(Auto_DBM)

    if Wind_type == 'Slow':
        ic(Wind_type)
        wind_array = np.clip(st.norm.rvs(
            370.530847, 88.585045, size=N), 1, 1000)
        gamma_array = np.clip(st.lognorm.rvs(
            0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08, size=N), 1.0e-08, 3.0e-7)

    else:
        ic(Wind_type)
        wind_array = np.clip(st.norm.rvs(
            579.057905, 67.870776, size=N), 1, 1000)
        gamma_array = np.clip(st.lognorm.rvs(
            0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08, size=N), 1.0e-08, 3.0e-7)

elif (PDBM == True) and (Auto_DBM == False):
    ic(PDBM)
    ic(Auto_DBM)

    print("Manual Entry of P-DBM Parameters are used for calculation.")
    print("Drag Parameter PDF is assumed to be Normal PDF")
    W = config['Solar_wind_Speed']
    dW = config['Error_Solar_wind_Speed']
    gamma = config['Drag_Parameter']
    dgamma = config['Error_Drag_Parameter']

    wind_array = np.random.normal(W, dW, N)
    gamma_array = np.clip(np.random.normal(gamma, dgamma, N), 1.0e-15, 3.0e-2)

elif (PDBM == False) and (Auto_DBM == True):
    ic(PDBM)
    ic(Auto_DBM)

    if Wind_type == 'Slow':
        ic(Wind_type)
        W = st.norm.mean(
            370.530847, 88.585045)
        gamma = st.lognorm.median(
            0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08)

    else:
        ic(Wind_type)
        W = st.norm.mean(
            579.057905, 67.870776)
        gamma = st.lognorm.median(
            0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08)

else:
    ic(PDBM)
    ic(Auto_DBM)
    print("Manual Entry of DBM Parameters are used for calculation.")
    W = config['Solar_wind_Speed']
    gamma = config['Drag_Parameter']


# =============================================================================
# Forecast Results for 1D (P)-DBM
# =============================================================================

if TwoD == False:
    ic(TwoD)

    if PDBM:
        T_array, V_array = f.PDBM(
            R0, dR0, R1, V0, dV0, gamma_array, wind_array, dT0, N)
        T_mean = np.nanmean(T_array)
        T_std = np.nanstd(T_array)
        V_mean = np.nanmean(V_array)
        V_std = np.nanstd(V_array)

        t_arrival = T0+(np.nanmedian(T_array)*3600)
        t_arrival_UTC = datetime.fromtimestamp(
            t_arrival).strftime("%Y-%m-%d %H:%M")
        T_arrival_UTC = datetime.strptime(t_arrival_UTC, "%Y-%m-%d %H:%M")

        print("1D P-DBM forecasting result.")
        print('------------------')
        print(f"* CME arrival forecast for {Target}.")
        print(
            f"* CME will arrive at {R_target:.3f} AU on(date and time) {T_arrival_UTC}.")
        print(
            f"* Mean transit time of CME is {T_mean:.3f} hrs and uncertainity is {T_std:.3f} hrs.")
        print(
            f"* Mean impact speed of CME is {V_mean:.3f} km/s and uncertainity is {V_std:.3f} km/s.")
        print('------------------')

    else:
        TT, VT = f.DBM(R0, R1, V0, gamma, W)

        t_arrival = T0+(TT*3600)
        t_arrival_UTC = datetime.fromtimestamp(
            t_arrival).strftime("%Y-%m-%d %H:%M")
        T_arrival_UTC = datetime.strptime(t_arrival_UTC, "%Y-%m-%d %H:%M")

        print("1D DBM forecasting result")
        print('------------------')
        print(f"* CME arrival forecast for {Target}.")
        print(
            f"* CME will arrive at {R_target:.3f} AU on(date and time) {T_arrival_UTC}")
        print(
            f"* Transit time of CME is {TT:.3f} hrs.")
        print(
            f"* Impact speed of CME is {VT:.3f} km/s.")
        print('------------------')
