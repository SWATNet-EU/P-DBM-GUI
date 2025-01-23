#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:10:50 2025

@author: astronish
"""

from PyQt6 import QtCore
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QMessageBox

from datetime import datetime
import functions as f
import astropy.units as u
import numpy as np
import scipy.stats as st
import ephem

from icecream import ic
import gc


class Actions2D:
    def __init__(self, ui):
        self.ui = ui
        self.setup_connections()

    def setup_connections(self):
        # Connecting Reset Buttons
        self.ui.Reset_2.clicked.connect(self.reset_2D)
        # Connect radio buttons to joint_conditions_2D
        self.ui.PDBM_Yes_2.toggled.connect(self.joint_conditions_2D)
        self.ui.PDBM_No_2.toggled.connect(self.joint_conditions_2D)
        self.ui.Auto_DBM_Yes_2.toggled.connect(self.joint_conditions_2D)
        self.ui.Auto_DBM_No_2.toggled.connect(self.joint_conditions_2D)
        self.ui.Calculate_2.clicked.connect(self.input_2D)

    def reset_2D(self):
        self.ui.CME_date_2.setDate(QtCore.QDate(2000, 1, 2))
        self.ui.CME_time_2.setTime(QtCore.QTime(0, 0, 0))
        self.reset_spin_boxes_2D()
        self.reset_combo_boxes_2D()
        self.reset_radio_buttons_2D()
        self.clear_results_2D()
        self.re_enable_widgets_2D()

    def reset_spin_boxes_2D(self):
        self.ui.CME_dt0_2.setValue(60.0)
        self.ui.R0_2.setValue(20.0)
        self.ui.dR0_2.setValue(1.0)
        self.ui.V0_2.setValue(1000.0)
        self.ui.dV0_2.setValue(100.0)
        self.ui.Omega.setValue(30.0)
        self.ui.dOmega.setValue(5.0)
        self.ui.Phi_CME.setValue(60.0)
        self.ui.dPhi_CME.setValue(5.0)
        self.ui.w_in_2.setValue(400.0)
        self.ui.dW_2.setValue(50.0)
        self.ui.Gamma_2.setValue(0.2)
        self.ui.dGamma_2.setValue(0.001)

    def reset_combo_boxes_2D(self):
        self.ui.Target_2.setCurrentIndex(2)  # Default to "Earth"
        self.ui.Cone_type.setCurrentIndex(1)  # Default to "Ice-Cream Cone"
        # Default to "Self-Similar Expansion"
        self.ui.Kinematic_Type.setCurrentIndex(0)

    def reset_radio_buttons_2D(self):
        self.ui.PDBM_Yes_2.setChecked(True)
        self.ui.Auto_DBM_Yes_2.setChecked(True)
        self.ui.W_slow_2.setChecked(True)
        self.ui.W_fast_2.setChecked(False)

    def clear_results_2D(self):
        self.ui.DBM_Result_2D.setText("")
        self.ui.RVT_plot_2.clear()
        self.ui.TT_Plot_2.clear()
        self.ui.V_plot_2.clear()
        self.ui.Plot.clear()

    def re_enable_widgets_2D(self):
        self.ui.dR0_2.setEnabled(True)
        self.ui.dV0_2.setEnabled(True)
        self.ui.CME_dt0_2.setEnabled(True)
        self.ui.W_slow_2.setEnabled(True)
        self.ui.W_fast_2.setEnabled(True)
        self.ui.w_in_2.setEnabled(True)
        self.ui.dW_2.setEnabled(True)
        self.ui.Gamma_2.setEnabled(True)
        self.ui.dGamma_2.setEnabled(True)
        self.ui.dPhi_CME.setEnabled(True)
        self.ui.dOmega.setEnabled(True)

    # Funtions to show avilable input based on the condition selected by user

    def joint_conditions_2D(self):
        if self.ui.PDBM_Yes_2.isChecked() and self.ui.Auto_DBM_Yes_2.isChecked():
            self.ui.dR0_2.setEnabled(True)
            self.ui.dV0_2.setEnabled(True)
            self.ui.CME_dt0_2.setEnabled(True)
            self.ui.dOmega.setEnabled(True)
            self.ui.dPhi_CME.setEnabled(True)
            self.ui.W_slow_2.setEnabled(True)
            self.ui.W_fast_2.setEnabled(True)
            self.ui.w_in_2.setEnabled(False)
            self.ui.dW_2.setEnabled(False)
            self.ui.Gamma_2.setEnabled(False)
            self.ui.dGamma_2.setEnabled(False)
        elif self.ui.PDBM_No_2.isChecked() and self.ui.Auto_DBM_Yes_2.isChecked():
            self.ui.dR0_2.setEnabled(False)
            self.ui.dV0_2.setEnabled(False)
            self.ui.CME_dt0_2.setEnabled(False)
            self.ui.dOmega.setEnabled(False)
            self.ui.dPhi_CME.setEnabled(False)
            self.ui.W_slow_2.setEnabled(True)
            self.ui.W_fast_2.setEnabled(True)
            self.ui.dW_2.setEnabled(False)
            self.ui.dGamma_2.setEnabled(False)
            self.ui.w_in_2.setEnabled(False)
            self.ui.Gamma_2.setEnabled(False)
        elif self.ui.PDBM_Yes_2.isChecked() and self.ui.Auto_DBM_No_2.isChecked():
            self.ui.dR0_2.setEnabled(True)
            self.ui.dV0_2.setEnabled(True)
            self.ui.CME_dt0_2.setEnabled(True)
            self.ui.dOmega.setEnabled(True)
            self.ui.dPhi_CME.setEnabled(True)
            self.ui.w_in_2.setEnabled(True)
            self.ui.Gamma_2.setEnabled(True)
            self.ui.dW_2.setEnabled(True)
            self.ui.dGamma_2.setEnabled(True)
            self.ui.W_slow_2.setEnabled(False)
            self.ui.W_fast_2.setEnabled(False)
        elif self.ui.PDBM_No_2.isChecked() and self.ui.Auto_DBM_No_2.isChecked():
            self.ui.dR0_2.setEnabled(False)
            self.ui.dV0_2.setEnabled(False)
            self.ui.CME_dt0_2.setEnabled(False)
            self.ui.dOmega.setEnabled(False)
            self.ui.dPhi_CME.setEnabled(False)
            self.ui.dW_2.setEnabled(False)
            self.ui.dGamma_2.setEnabled(False)
            self.ui.W_slow_2.setEnabled(False)
            self.ui.W_fast_2.setEnabled(False)
            self.ui.w_in_2.setEnabled(True)
            self.ui.Gamma_2.setEnabled(True)

    # 2D DBM input function

    def input_2D(self):
        """
        Clear variables from memory and then perfrom calculations.
        Below is cleaning part of the code.
        """
        
        # Spacific variables to clear
        variables_to_clear=["pdbm_choice", "auto_choice", "wind_type", "cme_date", "cme_time", "time_utc",
                            "T0","cme_dt0","r0","dr0", "v0", "dv0", "w_in", "dw", "gamma", "dgamma", "target",
                            "R_target", "r1", "wind_array", "gamma_array", "T_array", "V_array", "V0_array",
                            "plot", "T_PDF_plot", "V_PDF_plot"]
        
        for var in variables_to_clear:
            if hasattr(self, var):  # Check if the variable exists
               delattr(self, var)  # Delete the variable


        gc.collect()
        
       
        # Get values from radio buttons
        pdbm_choice = "Yes" if self.ui.PDBM_Yes_2.isChecked() else "No"
        auto_choice = "Yes" if self.ui.Auto_DBM_Yes_2.isChecked() else "No"
        wind_type = "Slow" if self.ui.W_slow_2.isChecked() else "Fast"

        # Get values from date and time inputs
        cme_date = self.ui.CME_date_2.date().toString("dd-MM-yyyy")
        cme_time = self.ui.CME_time_2.time().toString("hh:mm AP")

        # Combine date and time into a single string
        time_utc = f"{cme_date} {cme_time}"  # Example: "2000-01-01 00:00"

        # Convert the combined string into a datetime object
        # This variable convert the date time object to comatible format for the code.
        # converting string to date time object
        time_utc = datetime.strptime(time_utc, "%d-%m-%Y %I:%M %p")
        T0 = time_utc.timestamp()

        # Get values from spin boxes
        cme_dt0 = self.ui.CME_dt0_2.value()
        r0 = self.ui.R0_2.value()
        dr0 = self.ui.dR0_2.value()
        v0 = self.ui.V0_2.value()
        dv0 = self.ui.dV0_2.value()
        w_in = self.ui.w_in_2.value()
        dw = self.ui.dW_2.value()
        gamma = self.ui.Gamma_2.value()
        dgamma = self.ui.dGamma_2.value()

        # CME Angular Width and Propagation Direction
        omega = self.ui.Omega.value()
        domega = self.ui.dOmega.value()
        phi_cme = self.ui.Phi_CME.value()
        phi_cme = f.Phi_Correction(phi_cme, time_utc)
        dphi_cme = self.ui.dPhi_CME.value()

        # Get value from combo box and updating for calculation
        target = self.ui.Target_2.currentText()
        Phi_target, R_target = f.position(target, time_utc)

        # Check if R_target is NaN
        if np.isnan(R_target):
            # Show a popup window with an error message
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText(
                "Target is not there in the heliosphere during the provided time. \n All the input will be reset after you click on OK Button. ")
            msg.setWindowTitle("Invalid Input")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)

            # Execute the popup and handle the "OK" button
            if msg.exec() == QMessageBox.StandardButton.Ok:
                self.reset_2D()  # Reset inputs if "OK" is clicked
            return  # Exit the function to prevent further calculations
        r1 = ((R_target*u.au).to(u.km)).value

        cone_type = self.ui.Cone_type.currentText()
        kinematic_type = self.ui.Kinematic_Type.currentText()

        if cone_type == "Ice-Cream Cone" and omega > 90:
            # this is a limitating case of ice cream cone so we change the cone_type option for input
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText(
                "CME cone type with provided angular width is not possible.\n Hence, Concentric cone is considered for calculation.")
            msg.setWindowTitle("Invalid Input")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)

            # Execute the popup and handle the "OK" button
            if msg.exec() == QMessageBox.StandardButton.Ok:
                # Update the value in the variable
                cone_type = "Concentric Cone"
                # Update the value in the GUI
                self.ui.Cone_type.setCurrentText(cone_type)

        if cone_type == "Tangential Cone" and omega >= 90:
            # this is a limitating case of ice cream cone so we change the cone_type option for input
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText(
                "CME cone type with provided angular width is not possible.\n Hence, Concentric cone is considered for calculation.")
            msg.setWindowTitle("Invalid Input")
            msg.setStandardButtons(QMessageBox.StandardButton.Ok)

            # Execute the popup and handle the "OK" button
            if msg.exec() == QMessageBox.StandardButton.Ok:
                # Update the value in the variable
                cone_type = "Concentric Cone"
                # Update the value in the GUI
                self.ui.Cone_type.setCurrentText(cone_type)

        # Doing DBM calculation for 1D input
        if pdbm_choice == "Yes" and auto_choice == "Yes":

            # Number of Ensembles
            N = 5000

            if wind_type == 'Slow':
                wind_array = np.clip(st.norm.rvs(
                    370.530847, 88.585045, size=N), 1, 1000)
                gamma_array = np.clip(st.lognorm.rvs(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08, size=N), 1.0e-08, 3.0e-7)

            else:
                wind_array = np.clip(st.norm.rvs(
                    579.057905, 67.870776, size=N), 1, 1000)
                gamma_array = np.clip(st.lognorm.rvs(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08, size=N), 1.0e-08, 3.0e-7)


        elif pdbm_choice == "Yes" and auto_choice == "No":
            # Number of Ensembles
            N = 5000

            wind_array = np.random.normal(w_in, dw, N)
            gamma_array = np.clip(np.random.normal(
                gamma*1.0e-7, dgamma*1.0e-7, N), 1.0e-15, 3.0e-2)

            

        elif pdbm_choice == "No" and auto_choice == "Yes":
            if wind_type == 'Slow':
                W_2D = st.norm.mean(
                    370.530847, 88.585045)
                gamma_2D = st.lognorm.median(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08)
            else:
                W_2D = st.norm.mean(
                    579.057905, 67.870776)
                gamma_2D = st.lognorm.median(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08)
                
            if (phi_cme - omega) < Phi_target < (phi_cme + omega):
                alpha = np.abs(Phi_target - phi_cme)

                if cone_type == "Ice-Cream Cone" and kinematic_type == "Self-Similar Expansion":
                    TT, VT = f.DBM_IC_SSE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_IC_SSE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                    
                elif cone_type == "Ice-Cream Cone" and kinematic_type == "Flattening Cone Evolution":
                    TT, VT = f.DBM_IC_FCE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_IC_FCE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                    
                elif cone_type == "Tangential Cone" and kinematic_type == "Flattening Cone Evolution":
                    TT, VT = f.DBM_TC_FCE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_TC_FCE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                    
                elif cone_type == "Tangential Cone" and kinematic_type == "Self-Similar Expansion":
                    TT, VT = f.DBM_TC_SSE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_TC_SSE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                elif cone_type == "Concentric Cone": # same as 1D DBM
                    TT, VT = f.DBM(r0, r1, v0, gamma_2D, W_2D)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_RVT_plot(time_utc, TT*3600, r0,
                                          v0, gamma_2D, W_2D, R_target, tdate)

                output_text = (f"* CME hits the target.<br>"
                               f"* Space Weather Alert. <br>"
                               f"* CME arrival forecast for <b>{target}</b>.<br>"
                               f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                               f"* Mean transit time of CME is <b>{TT:.3f} hrs</b>.<br>"
                               f"* Mean impact speed of CME is <b>{VT:.3f} km/s</b>.<br>"
                               f"* Median value of drag parameter γ is <b>{gamma_2D/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                               f"* Median Value of solar wind speed w is <b>{W_2D:.2f} km/s</b>.</br> ")

            else:
                TT, VT = f.DBM(r0, r1, v0, gamma_2D, W_2D)

                t_arrival = T0+(TT*3600)
                t_arrival_UTC = datetime.fromtimestamp(
                    t_arrival).strftime("%Y-%m-%d %H:%M")
                tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                T_arrival_UTC = datetime.strptime(
                    t_arrival_UTC, "%Y-%m-%d %H:%M")

                # Making Plots
                plot = f.DBM_RVT_plot(time_utc, TT*3600, r0,
                                      v0, gamma_2D, W_2D, R_target, tdate)

                output_text = (f"* Yay !!! CME misses the target.<br>"
                               f"* Model Calculation for Research Purpose.<br>"
                               f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                               f"* Mean transit time of CME is <b>{TT:.3f} hrs</b>.<br>"
                               f"* Mean impact speed of CME is <b>{VT:.3f} km/s</b>.<br>"
                               f"* Median value of drag parameter γ is <b>{gamma_2D/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                               f"* Median Value of solar wind speed w is <b>{W_2D:.2f} km/s</b>.</br> ")

            # Showing Output
            self.ui.DBM_Result_2D.show()
            self.ui.DBM_Result_2D_scroll.show()  # Show the scroll area
            self.ui.DBM_Result_2D_scroll.setGeometry(390, 60, 400, 300)
            self.ui.DBM_Result_2D.clear()
            self.ui.DBM_Result_2D.setText(output_text)
            # self.ui.DBM_Result_2D.adjustSize()
            # self.ui.DBM_Result_2D.setFixedWidth(501)
            # self.ui.DBM_Result_2D.setWordWrap(True)  # Enable word wrapping
            # Showing Plot
            self.ui.TT_Plot_2.hide()
            self.ui.V_plot_2.hide()
            pixmap = QPixmap()
            pixmap.loadFromData(plot.read())
            self.ui.RVT_plot_2.show()
            self.ui.RVT_plot_2.clear()
            self.ui.RVT_plot_2.setGeometry(390, 400, 600, 400)
            self.ui.RVT_plot_2.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.RVT_plot_2.setScaledContents(True)

            

        elif pdbm_choice == "No" and auto_choice == "No":
            W_2D = w_in
            gamma_2D = gamma*1.0e-7

            if (phi_cme - omega) < Phi_target < (phi_cme + omega):
                alpha = np.abs(Phi_target - phi_cme)

                if cone_type == "Ice-Cream Cone" and kinematic_type == "Self-Similar Expansion":
                    TT, VT = f.DBM_IC_SSE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)
                    ic(cone_type)
                    ic(kinematic_type)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_IC_SSE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                elif cone_type == "Ice-Cream Cone" and kinematic_type == "Flattening Cone Evolution":
                    TT, VT = f.DBM_IC_FCE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)
                    ic(cone_type)
                    ic(kinematic_type)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_IC_FCE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                elif cone_type == "Tangential Cone" and kinematic_type == "Flattening Cone Evolution":
                    TT, VT = f.DBM_TC_FCE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)
                    ic(cone_type)
                    ic(kinematic_type)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_TC_FCE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                elif cone_type == "Tangential Cone" and kinematic_type == "Self-Similar Expansion":
                    TT, VT = f.DBM_TC_SSE(
                        r0, r1, v0, gamma_2D, W_2D, omega, alpha)
                    ic(cone_type)
                    ic(kinematic_type)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_2D_RVT_TC_SSE_plot(time_utc, TT*3600, r0,
                                                    v0, gamma_2D, W_2D, R_target, tdate, omega, alpha)
                elif cone_type == "Concentric Cone": # same as 1D DBM
                    TT, VT = f.DBM(r0, r1, v0, gamma_2D, W_2D)
                    ic(cone_type)
                    ic(kinematic_type)

                    t_arrival = T0+(TT*3600)
                    t_arrival_UTC = datetime.fromtimestamp(
                        t_arrival).strftime("%Y-%m-%d %H:%M")
                    tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                    T_arrival_UTC = datetime.strptime(
                        t_arrival_UTC, "%Y-%m-%d %H:%M")
                    # Making Plots
                    plot = f.DBM_RVT_plot(time_utc, TT*3600, r0,
                                          v0, gamma_2D, W_2D, R_target, tdate)

                output_text = (f"* CME hits the target.<br>"
                               f"* Space Weather Alert. <br>"
                               f"* CME arrival forecast for <b>{target}</b>.<br>"
                               f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                               f"* Mean transit time of CME is <b>{TT:.3f} hrs</b>.<br>"
                               f"* Mean impact speed of CME is <b>{VT:.3f} km/s</b>.<br>"
                               f"* Median value of drag parameter γ is <b>{gamma_2D/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                               f"* Median Value of solar wind speed w is <b>{W_2D:.2f} km/s</b>.</br> ")

            else:
                TT, VT = f.DBM(r0, r1, v0, gamma_2D, W_2D)
                ic()

                t_arrival = T0+(TT*3600)
                t_arrival_UTC = datetime.fromtimestamp(
                    t_arrival).strftime("%Y-%m-%d %H:%M")
                tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
                T_arrival_UTC = datetime.strptime(
                    t_arrival_UTC, "%Y-%m-%d %H:%M")

                # Making Plots
                plot = f.DBM_RVT_plot(time_utc, TT*3600, r0,
                                      v0, gamma_2D, W_2D, R_target, tdate)

                output_text = (f"* Yay !!! CME misses the target.<br>"
                               f"* Model Calculation for Research Purpose.<br>"
                               f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                               f"* Mean transit time of CME is <b>{TT:.3f} hrs</b>.<br>"
                               f"* Mean impact speed of CME is <b>{VT:.3f} km/s</b>.<br>"
                               f"* Median value of drag parameter γ is <b>{gamma_2D/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                               f"* Median Value of solar wind speed w is <b>{W_2D:.2f} km/s</b>.</br> ")

            # Showing Output
            self.ui.DBM_Result_2D.show()
            self.ui.DBM_Result_2D_scroll.show()  # Show the scroll area
            self.ui.DBM_Result_2D_scroll.setGeometry(390, 60, 400, 300)
            self.ui.DBM_Result_2D.clear()
            self.ui.DBM_Result_2D.setText(output_text)
            # self.ui.DBM_Result_2D.adjustSize()
            # self.ui.DBM_Result_2D.setFixedWidth(501)
            # self.ui.DBM_Result_2D.setWordWrap(True)  # Enable word wrapping
            # Showing Plot
            self.ui.TT_Plot_2.hide()
            self.ui.V_plot_2.hide()
            pixmap = QPixmap()
            pixmap.loadFromData(plot.read())
            self.ui.RVT_plot_2.show()
            self.ui.RVT_plot_2.clear()
            self.ui.RVT_plot_2.setGeometry(390, 400, 600, 400)
            self.ui.RVT_plot_2.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.RVT_plot_2.setScaledContents(True)
