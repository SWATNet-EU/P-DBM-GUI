#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 15:33:22 2024

@author: astronish16
"""


# Numeric imports
import sys
from icecream import ic
import numpy as np
import scipy as sc
from scipy.optimize import newton
from astroquery.jplhorizons import Horizons
from astropy.time import Time
import astropy.units as u
import ephem
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime, date, timedelta
from io import BytesIO


# Some Lists and Array for automization purposes
Inner_Planets = ["Mercury", "Venus", "Earth", "Mars"]
Outer_Planets = ["Jupiter", "Saturn", "Uranus", "Neptune"]
Space_Crafts = ["Messenger", "VEX", "PSP", "SolO", "BepiCol", "Spitzer",
                "Wind", "ST-A", "ST-B", "Kepler", "Ulysses", "MSL", "Maven", "Juno"]


objects_list = Inner_Planets + Outer_Planets + Space_Crafts


# Function required to find DBM solution
def func(t, r0, r1, v0, gamma, w):

    if v0 >= w:
        gamma = gamma
    else:
        # only possible contion is v0<w.
        gamma = -1*gamma

    p1 = 1 + gamma*(v0-w)*t
    y = -r1 + r0 + w*t + (np.log(p1)/gamma)
    y1 = w + ((v0-w)/p1)

    return y, y1


# Wrapper function for y and y1.
'''
It is necessary for the DBM function.
If other solution is there then it need to be find.
calleable function can enhance the performance
'''


def func_y(t, r0, r1, v0, gamma, w):
    y, _ = func(t, r0, r1, v0, gamma, w)
    return y


def func_y1(t, r0, r1, v0, gamma, w):
    _, y1 = func(t, r0, r1, v0, gamma, w)
    return y1

# Function to plot RVT plot for DBM


def RV(t, r0, v0, gamma, w):

    if v0 >= w:
        gamma = gamma
    else:
        # only possible contion is v0<w.
        gamma = -1*gamma

    p1 = 1 + gamma*(v0-w)*t
    r = r0 + w*t + (np.log(p1)/gamma)
    v = w + ((v0-w)/p1)

    return r, v


def DBM_RVT_plot(time_utc, TT, r0, v0, gamma, w, r_target, tdate):

    dt = 3600     # unit is second
    t_ary = np.arange(0, TT*1.1*3600, 80)
    Time = [time_utc + timedelta(seconds=i) for i in t_ary]

    R, V = RV(t_ary, r0, v0, gamma, w)
    R = (R*u.km).to(u.R_sun).value
    r_target = (r_target*u.au).to(u.R_sun).value

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.grid()
    color = 'tab:red'
    ax1.set_xlabel('time (UTC date hour)', fontsize=17)
    ax1.set_ylabel('R (solar radius)', color=color, fontsize=17)
    ax1.plot(Time, R, color=color, label="Distance")
    ax1.axvline(tdate, linestyle='--', color="black", label="Arrival Time")
    ax1.axhline(r_target, label=f"R_target = {r_target:.2f}")
    ax1.tick_params(axis='y', labelcolor=color, labelsize=17)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:cyan'
    # we already handled the x-label with ax1
    ax2.set_ylabel('V (km/s)', color=color, fontsize=17)
    ax2.plot(Time, V, color=color, label="Speed")
    ax2.tick_params(axis='y', labelcolor=color,labelsize=17)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()

    lines = lines_1 + lines_2
    labels = labels_1 + labels_2

    ax1.legend(lines, labels, fontsize=17, loc=3)
    plt.grid(True)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


def PDBM_RVT_plot(time_utc, TT_array, r0, v0_array, gamma_array, w_array, r_target, tdate):

    TT = np.nanmedian(TT_array)
    t_ary = np.arange(0, TT*1.1*3600, 80)
    Time = [time_utc + timedelta(seconds=i) for i in t_ary]

    i_array = np.arange(0, len(v0_array), 1)
    R_matrix = np.zeros((len(v0_array), len(t_ary)))
    V_matrix = np.zeros((len(v0_array), len(t_ary)))

    for v0, w, g, i in zip(v0_array, w_array, gamma_array, i_array):
        y, y1 = RV(t_ary, r0, v0, g, w)
        y = (y*u.km).to(u.R_sun).value
        R_matrix[i] = y
        V_matrix[i] = y1

    r_target = (r_target*u.au).to(u.R_sun).value
    R_median = np.nanmedian(R_matrix, axis=0)
    V_median = np.nanmedian(V_matrix, axis=0)
    R_max = np.max(R_matrix, axis=0)
    R_min = np.min(R_matrix, axis=0)
    V_max = np.max(V_matrix, axis=0)
    V_min = np.min(V_matrix, axis=0)

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.grid()
    color = 'tab:red'
    ax1.set_xlabel('time (UTC date hour)', fontsize=17)
    ax1.set_ylabel('R (solar radius)', color=color, fontsize=17)
    ax1.plot(Time, R_median, color=color, label="Distance")
    ax1.fill_between(Time, R_max, R_min, alpha=0.25,
                     linewidth=0, color=color)
    ax1.axvline(tdate, linestyle='--', color="black", label="Arrival Time")
    ax1.axhline(r_target, label=f"R_target = {r_target:.2f}")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:cyan'
    # we already handled the x-label with ax1
    ax2.set_ylabel('V (km/s)', color=color, fontsize=17)
    ax2.plot(Time, V_median, color=color, label="Speed")
    ax2.fill_between(Time, V_max, V_min, alpha=0.25,
                     linewidth=0, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()

    lines = lines_1 + lines_2
    labels = labels_1 + labels_2

    ax1.legend(lines, labels, fontsize=17, loc=3)
    plt.grid(True)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


def TT_plot(T):
    mean = np.nanmean(T)
    median = np.nanmedian(T)
    plt.style.use('seaborn-v0_8')
    plt.hist(T, bins=50, density=True)
    plt.xlim(np.nanmin(T), np.nanmax(T))
    plt.axvline(mean, color='red', label=f"Mean: {mean:.2f} hr ")
    plt.axvline(median, color='black', label=f"Median: {median:.2f} hr ")
    plt.title("Transit Time Distribution")
    plt.xlabel("Transit Time (hrs)")
    plt.xlim(0,np.nanmax(T))
    plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


def V_plot(V):
    mean = np.nanmean(V)
    median = np.nanmedian(V)
    plt.style.use('seaborn-v0_8')
    plt.hist(V, bins=50, density=True)
    plt.xlim(np.nanmin(V), np.nanmax(V))
    plt.axvline(mean, color='red', label=f"Mean: {mean:.2f} km/s ")
    plt.axvline(median, color='black', label=f"Median: {median:.2f} km/s ")
    plt.title("Arrival Speed Distribution")
    plt.xlabel("Arrival Speed (km/s)")
    plt.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


# Function to find DBM solution.
'''
This function provides transist time and impact speed
of CME under DBM approximation
'''


def DBM(r0, r1, v0, gamma, w):
    t1 = newton(func=func_y, fprime=func_y1, x0=30*3600,
                args=[r0, r1, v0, gamma, w], disp=False, maxiter=30)

    dv = v0-w
    p1 = 1 + (gamma*np.abs(dv)*t1)
    v1 = w + (dv/p1)

    # Transit time is in hours and impact speed is in km/s.
    return t1/3600.0, v1


def PDBM(r0, dr0, r1, v0, dv0, gamma_array, wind_array, dt0, N):
    r0_array = np.random.normal(r0, dr0/3.0, N)
    v0_array = np.random.normal(v0, dv0/3.0, N)
    r1_array = np.random.normal(r1, 0.05*r1/3.0, N) # including 5% error in the target distace
    t0_array = np.random.normal(0, dt0/(60.0*3), N)

    # Arrays to store Output
    TT_array = np.zeros_like(r0_array)
    V_array = np.zeros_like(r0_array)

    for i in range(0, N):
        TT_array[i], V_array[i] = DBM(
            r0_array[i], r1_array[i], v0_array[i], gamma_array[i], wind_array[i])
        TT_array[i] = TT_array[i]+t0_array[i]

    return TT_array, V_array, v0_array


# Function for 2D DBM
'''
while moving to 2D version of model,one has to consider two important variable.
(1) CME Propagation Direction
(2) Target longitude.
Difference of these two quantity is called alpha.
Also, CME propagation direction has been measured withrespect to the Earth in cone geometry.
Therefore, we also need to consider coordinate system.
'''

# Correction in central meridian as per Heliocentric Ecliptic coordinate system.
'''
We are doing this because we are using JPL Horiozn for ephemeris to determine target information.
In JPL Horizon Heliocentric ecliptic coordinate system is used so position of Earth is not Fixed.
While in cone model Heliocentric Stonyhurst coordinate system is used where Sun-Earth line is always
correspond to the 0$^o$ longitude.
'''


def Phi_Correction(phi_cme, time_utc):
    earth = ephem.Sun()
    earth.compute(time_utc)
    phi_corrected = np.rad2deg(np.deg2rad(phi_cme) + earth.hlon)
    return phi_corrected


'''
For 2D cone, There are 3 possiblities (see Schwenn et al, 2005)
(1) ICME leading edge is concentric arc with solar surface.
    The application of this geometry is same as 1D DBM.
(2) ICME leading edge is semi circle.
    This geometry looks like a ice cream cone
(3) ICME leading edge is circular arc and tangentially connect to the ICME legs.
    Application of this geometry is bit difficult.
'''

'''
When we consider a geometry, possiblity of two different type of evolution is arise.
(1) Self Similar Expansion: CME maintain it's shape during propagation.
(2) Flattening Cone Evolution: Each and every point on CME edge follows DBM.
For more detailed informatio: Check the documantation.
'''

# Function to calculate speed and distance at alpha angle

# ICE Cream Cone only


def IC_RV_alpha(omega, alpha, r0, v0):
    omega = np.deg2rad(omega)
    alpha = np.deg2rad(alpha)
    r01 = r0 * (np.cos(alpha) + ((np.tan(omega))**2 -
                (np.sin(alpha))**2)**0.5)/(1 + np.tan(omega))
    v01 = v0 * (np.cos(alpha) + ((np.tan(omega))**2 -
                (np.sin(alpha))**2)**0.5)/(1 + np.tan(omega))
    return r01, v01


def IC_R_alpha_inv(omega, alpha, r1):
    omega = np.deg2rad(omega)
    alpha = np.deg2rad(alpha)
    r1_apex = r1 * (1 + np.tan(omega))/((np.cos(alpha)) +
                                        (((np.tan(omega))**2.0 - (np.sin(alpha))**2.0)**0.5))
    return r1_apex


def TC_RV_alpha(omega, alpha, r0, v0):
    omega = np.deg2rad(omega)
    alpha = np.deg2rad(alpha)
    r01 = r0 * (np.cos(alpha) + ((np.sin(omega))**2 -
                (np.sin(alpha))**2)**0.5)/(1 + np.sin(omega))
    v01 = v0 * (np.cos(alpha) + ((np.sin(omega))**2 -
                (np.sin(alpha))**2)**0.5)/(1 + np.sin(omega))
    return r01, v01


def TC_R_alpha_inv(omega, alpha, r1):
    omega = np.deg2rad(omega)
    alpha = np.deg2rad(alpha)
    r1_apex = r1 * (1 + np.sin(omega))/((np.cos(alpha)) +
                                        (((np.sin(omega))**2.0 - (np.sin(alpha))**2.0)**0.5))
    return r1_apex


def DBM_2D_RVT_IC_SSE_plot(time_utc, TT, r0, v0, gamma, w, r_target, tdate, omega, alpha):

    t_ary = np.arange(0, TT*1.1*3600, 80)
    Time = [time_utc + timedelta(seconds=i) for i in t_ary]

    R, V = RV(t_ary, r0, v0, gamma, w)
    R = (R*u.km).to(u.R_sun).value
    r_target = (r_target*u.au).to(u.R_sun).value

    R_ary, V_ary = IC_RV_alpha(omega, alpha, R, V)

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.grid()
    color = 'tab:red'
    ax1.set_xlabel('time (UTC date hour)', fontsize=17)
    ax1.set_ylabel('R (solar radius)', color=color, fontsize=17)
    ax1.plot(Time, R_ary, color=color, label="Distance")
    ax1.axvline(tdate, linestyle='--', color="black", label="Arrival Time")
    ax1.axhline(r_target, label=f"R_target = {r_target:.2f}")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:cyan'
    # we already handled the x-label with ax1
    ax2.set_ylabel('V (km/s)', color=color, fontsize=17)
    ax2.plot(Time, V_ary, color=color, label="Speed")
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()

    lines = lines_1 + lines_2
    labels = labels_1 + labels_2

    ax1.legend(lines, labels, fontsize=17, loc=3)
    plt.grid(True)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


def DBM_2D_RVT_IC_FCE_plot(time_utc, TT, r0, v0, gamma, w, r_target, tdate, omega, alpha):

    t_ary = np.arange(0, TT*1.1*3600, 80)
    Time = [time_utc + timedelta(seconds=i) for i in t_ary]

    R0_a, V0_a = IC_RV_alpha(omega, alpha, r0, v0)

    R, V = RV(t_ary, R0_a, V0_a, gamma, w)
    R = (R*u.km).to(u.R_sun).value
    r_target = (r_target*u.au).to(u.R_sun).value

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.grid()
    color = 'tab:red'
    ax1.set_xlabel('time (UTC date hour)', fontsize=17)
    ax1.set_ylabel('R (solar radius)', color=color, fontsize=17)
    ax1.plot(Time, R, color=color, label="Distance")
    ax1.axvline(tdate, linestyle='--', color="black", label="Arrival Time")
    ax1.axhline(r_target, label=f"R_target = {r_target:.2f}")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:cyan'
    # we already handled the x-label with ax1
    ax2.set_ylabel('V (km/s)', color=color, fontsize=17)
    ax2.plot(Time, V, color=color, label="Speed")
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()

    lines = lines_1 + lines_2
    labels = labels_1 + labels_2

    ax1.legend(lines, labels, fontsize=17, loc=3)
    plt.grid(True)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


def DBM_2D_RVT_TC_SSE_plot(time_utc, TT, r0, v0, gamma, w, r_target, tdate, omega, alpha):

    t_ary = np.arange(0, TT*1.1*3600, 80)
    Time = [time_utc + timedelta(seconds=i) for i in t_ary]

    R, V = RV(t_ary, r0, v0, gamma, w)
    R = (R*u.km).to(u.R_sun).value
    r_target = (r_target*u.au).to(u.R_sun).value

    R_ary, V_ary = TC_RV_alpha(omega, alpha, R, V)

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.grid()
    color = 'tab:red'
    ax1.set_xlabel('time (UTC date hour)', fontsize=17)
    ax1.set_ylabel('R (solar radius)', color=color, fontsize=17)
    ax1.plot(Time, R_ary, color=color, label="Distance")
    ax1.axvline(tdate, linestyle='--', color="black", label="Arrival Time")
    ax1.axhline(r_target, label=f"R_target = {r_target:.2f}")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:cyan'
    # we already handled the x-label with ax1
    ax2.set_ylabel('V (km/s)', color=color, fontsize=17)
    ax2.plot(Time, V_ary, color=color, label="Speed")
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()

    lines = lines_1 + lines_2
    labels = labels_1 + labels_2

    ax1.legend(lines, labels, fontsize=17, loc=3)
    plt.grid(True)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


def DBM_2D_RVT_TC_FCE_plot(time_utc, TT, r0, v0, gamma, w, r_target, tdate, omega, alpha):

    t_ary = np.arange(0, TT*1.1*3600, 80)
    Time = [time_utc + timedelta(seconds=i) for i in t_ary]

    R0_a, V0_a = TC_RV_alpha(omega, alpha, r0, v0)

    R, V = RV(t_ary, R0_a, V0_a, gamma, w)
    R = (R*u.km).to(u.R_sun).value
    r_target = (r_target*u.au).to(u.R_sun).value

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.grid()
    color = 'tab:red'
    ax1.set_xlabel('time (UTC date hour)', fontsize=17)
    ax1.set_ylabel('R (solar radius)', color=color, fontsize=17)
    ax1.plot(Time, R, color=color, label="Distance")
    ax1.axvline(tdate, linestyle='--', color="black", label="Arrival Time")
    ax1.axhline(r_target, label=f"R_target = {r_target:.2f}")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:cyan'
    # we already handled the x-label with ax1
    ax2.set_ylabel('V (km/s)', color=color, fontsize=17)
    ax2.plot(Time, V, color=color, label="Speed")
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()

    lines = lines_1 + lines_2
    labels = labels_1 + labels_2

    ax1.legend(lines, labels, fontsize=17, loc=3)
    plt.grid(True)

    # Save the plot to an in-memory buffer
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)  # Move to the start of the buffer
    plt.close()  # Close the plot to free resources

    return buffer


'''
in self similar expansion we need to find the position of apex first.
Then we need to solve the boundary condition accordingly,
as target distance could be anywhere along the ICME leading edge.
With this corrected boundary condition we get a Travel time and arrival speed of CME apex.
Travel time doesn't change with angle correction.
But arrival spped do chnage.
So, we need to find the cosine component of CME apex speed which arrival speed of CME at target.
'''


def DBM_IC_SSE(r0, r1, v0, gamma, w, omega, alpha):
    R1_apex = IC_R_alpha_inv(omega, alpha, r1)
    TT, V1_apex = DBM(r0, R1_apex, v0, gamma, w)
    _, V1 = IC_RV_alpha(omega, alpha, R1_apex, V1_apex)
    return TT, V1


def DBM_TC_SSE(r0, r1, v0, gamma, w, omega, alpha):
    R1_apex = TC_R_alpha_inv(omega, alpha, r1)
    TT, V1_apex = DBM(r0, R1_apex, v0, gamma, w)
    _, V1 = TC_RV_alpha(omega, alpha, R1_apex, V1_apex)
    return TT, V1


def DBM_IC_FCE(r0, r1, v0, gamma, w, omega, alpha):
    r0_a, v0_a = IC_RV_alpha(omega, alpha, r0, v0)
    TT, V1 = DBM(r0_a, r1, v0_a, gamma, w)
    return TT, V1


def DBM_TC_FCE(r0, r1, v0, gamma, w, omega, alpha):
    r0_a, v0_a = TC_RV_alpha(omega, alpha, r0, v0)
    TT, V1 = DBM(r0_a, r1, v0_a, gamma, w)
    return TT, V1


# def DBM_Self_Similar_Expansion(r0, r1, v0, gamma, w, omega, alpha):
#     R1_apex = R_alpha_inv(omega, alpha, r1)
#     TT, V1_apex = DBM(r0, R1_apex, v0, gamma, w)
#     _, V1 = RV_alpha(omega, alpha, R1_apex, V1_apex)
#     return TT, V1


# def Forecast_SSE(r0, r1, v0, gamma, w, omega, phi_cme, phi_target):
#     if (phi_cme - omega) <= phi_target <= (phi_cme + omega):
#         print("Ohh no! CME hits the target")
#         print("Space Weather Alert")
#         alpha = np.abs(phi_cme - phi_target)
#         TT, V1 = DBM_Self_Similar_Expansion(r0, r1, v0, gamma, w, omega, alpha)
#     else:
#         print("Yay !!! CME misses the target")
#         print("Model Calculation for Research Purpose")
#         print("1D-DBM Model values")
#         TT, V1 = DBM(r0, r1, v0, gamma, w)
#     return TT, V1


# def DBM_Flattening_Cone(r0, r1, v0, gamma, w, omega, alpha):
#     r0_a, v0_a = RV_alpha(omega, alpha, r0, v0)
#     TT, V1 = DBM(r0_a, r1, v0_a, gamma, w)
#     return TT, V1


# def Forecast_FC(r0, r1, v0, gamma, w, omega, phi_cme, phi_target):
#     if (phi_cme - omega) <= phi_target <= (phi_cme + omega):
#         print("Ohh no! CME hits the target")
#         print("Space Weather Alert")
#         alpha = np.abs(phi_cme - phi_target)
#         TT, V1 = DBM_Flattening_Cone(r0, r1, v0, gamma, w, omega, alpha)
#     else:
#         print("Yay !!! CME misses the target")
#         print("Model Calculation for Research Purpose")
#         print("1D-DBM Model values")
#         TT, V1 = DBM(r0, r1, v0, gamma, w)
#     return TT, V1


def horizon_id(target):
    if target not in objects_list:
        print("Provided heliospheric object doesn't exist in the object list")
        print("Please check the object list to find correct object name.")
        print(f"Available Objects: {objects_list}")
        sys.exit()
    elif target == 'Mercury':
        h_id = 199
    elif target == 'Venus':
        h_id = 299
    elif target == 'Earth':
        h_id = 399
    elif target == 'Mars':
        h_id = 499
    elif target == 'Jupiter':
        h_id = 599
    elif target == 'Saturn':
        h_id = 699
    elif target == 'Uranus':
        h_id = 799
    elif target == 'Neptune':
        h_id == 899
    elif target == 'Messenger':
        h_id = -236
    elif target == 'VEX':
        h_id = -248
    elif target == 'PSP':
        h_id = -96
    elif target == 'SolO':
        h_id = -144
    elif target == 'BepiCol':
        h_id = -121
    elif target == 'Spitzer':
        h_id = -79
    elif target == 'Wind':
        h_id = -8
    elif target == 'ST-A':
        h_id = -234
    elif target == 'ST-B':
        h_id = -235
    elif target == 'Kepler':
        h_id = -227
    elif target == 'Ulysses':
        h_id = -55
    elif target == 'MSL':
        h_id = -76
    elif target == 'Maven':
        h_id = -202
    elif target == 'Juno':
        h_id = -61

    return str(h_id)


def position(target, date):
    tjd = Time(date).jd
    h_id = horizon_id(target)
    try:
        obj = Horizons(id=h_id, location='@sun', epochs=tjd)
        obj_eph = obj.ephemerides()
        obj_name = obj_eph['targetname']
        phi_obj = obj_eph['EclLon'][0]
        r_obj = obj_eph['r'][0]
    except:
        print(f"No {target} during the provided date")
        phi_obj = np.nan
        r_obj = np.nan
        
    return phi_obj, r_obj
