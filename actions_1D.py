#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:54:15 2024

@author: smp21rhm
"""

from PyQt6 import QtCore
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import QMessageBox

from datetime import datetime
import functions as f
from astropy.constants import R_sun
import astropy.units as u
import numpy as np
import scipy.stats as st
import ephem

from icecream import ic
import gc


class Actions1D:
    def __init__(self, ui):
        self.ui = ui
        self.setup_connections()

    def setup_connections(self):
        # Connecting Reset Buttons
        self.ui.Reset.clicked.connect(self.reset_1D)
        # Connect radio buttons to joint_conditions_1D
        self.ui.PDBM_Yes.toggled.connect(self.joint_conditions_1D)
        self.ui.PDBM_No.toggled.connect(self.joint_conditions_1D)
        self.ui.Auto_DBM_Yes.toggled.connect(self.joint_conditions_1D)
        self.ui.Auto_DBM_No.toggled.connect(self.joint_conditions_1D)
        self.ui.Calculate.clicked.connect(self.input_1D)

    def reset_1D(self):
        self.ui.CME_date.setDate(QtCore.QDate(2000, 1, 2))
        self.ui.CME_time.setTime(QtCore.QTime(0, 0, 0))
        self.reset_spin_boxes_1D()
        self.reset_combo_boxes_1D()
        self.reset_radio_buttons_1D()
        self.clear_results_1D()
        self.re_enable_widgets_1D()

    def reset_spin_boxes_1D(self):
        self.ui.CME_dt0.setValue(60.0)
        self.ui.R0.setValue(20.0)
        self.ui.dR0.setValue(1.0)
        self.ui.V0.setValue(1000.0)
        self.ui.dV0.setValue(100.0)
        self.ui.w_in.setValue(400.0)
        self.ui.dW.setValue(50.0)
        self.ui.Gamma.setValue(0.2)
        self.ui.dGamma.setValue(0.001)

    def reset_combo_boxes_1D(self):
        self.ui.Target.setCurrentIndex(2)  # Default to "Earth"

    def reset_radio_buttons_1D(self):
        self.ui.PDBM_Yes.setChecked(True)
        self.ui.Auto_DBM_Yes.setChecked(True)
        self.ui.W_slow.setChecked(True)
        self.ui.W_fast.setChecked(False)

    def clear_results_1D(self):
        self.ui.DBM_Result_1D.setText("")
        self.ui.RVT_plot.clear()
        self.ui.TT_Plot.clear()
        self.ui.V_plot.clear()

    def re_enable_widgets_1D(self):
        self.ui.dR0.setEnabled(True)
        self.ui.dV0.setEnabled(True)
        self.ui.CME_dt0.setEnabled(True)
        self.ui.W_slow.setEnabled(True)
        self.ui.W_fast.setEnabled(True)
        self.ui.w_in.setEnabled(True)
        self.ui.dW.setEnabled(True)
        self.ui.Gamma.setEnabled(True)
        self.ui.dGamma.setEnabled(True)

    # Funtions to show avilable input based on the condition selected by user

    def joint_conditions_1D(self):
        if self.ui.PDBM_Yes.isChecked() and self.ui.Auto_DBM_Yes.isChecked():
            self.ui.dR0.setEnabled(True)
            self.ui.dV0.setEnabled(True)
            self.ui.CME_dt0.setEnabled(True)
            self.ui.W_slow.setEnabled(True)
            self.ui.W_fast.setEnabled(True)
            self.ui.w_in.setEnabled(False)
            self.ui.dW.setEnabled(False)
            self.ui.Gamma.setEnabled(False)
            self.ui.dGamma.setEnabled(False)
        elif self.ui.PDBM_No.isChecked() and self.ui.Auto_DBM_Yes.isChecked():
            self.ui.dR0.setEnabled(False)
            self.ui.dV0.setEnabled(False)
            self.ui.CME_dt0.setEnabled(False)
            self.ui.W_slow.setEnabled(True)
            self.ui.W_fast.setEnabled(True)
            self.ui.dW.setEnabled(False)
            self.ui.dGamma.setEnabled(False)
            self.ui.w_in.setEnabled(False)
            self.ui.Gamma.setEnabled(False)
        elif self.ui.PDBM_Yes.isChecked() and self.ui.Auto_DBM_No.isChecked():
            self.ui.dR0.setEnabled(True)
            self.ui.dV0.setEnabled(True)
            self.ui.CME_dt0.setEnabled(True)
            self.ui.w_in.setEnabled(True)
            self.ui.Gamma.setEnabled(True)
            self.ui.dW.setEnabled(True)
            self.ui.dGamma.setEnabled(True)
            self.ui.W_slow.setEnabled(False)
            self.ui.W_fast.setEnabled(False)
        elif self.ui.PDBM_No.isChecked() and self.ui.Auto_DBM_No.isChecked():
            self.ui.dR0.setEnabled(False)
            self.ui.dV0.setEnabled(False)
            self.ui.CME_dt0.setEnabled(False)
            self.ui.dW.setEnabled(False)
            self.ui.dGamma.setEnabled(False)
            self.ui.W_slow.setEnabled(False)
            self.ui.W_fast.setEnabled(False)
            self.ui.w_in.setEnabled(True)
            self.ui.Gamma.setEnabled(True)

    # 1D DBM input function

    def input_1D(self):
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
        
        # Retrieve input values from 1D-DBM Tab

        # Get values from radio buttons
        pdbm_choice = "Yes" if self.ui.PDBM_Yes.isChecked() else "No"
        auto_choice = "Yes" if self.ui.Auto_DBM_Yes.isChecked() else "No"
        wind_type = "Slow" if self.ui.W_slow.isChecked() else "Fast"

        # Get values from date and time inputs
        cme_date = self.ui.CME_date.date().toString("dd-MM-yyyy")
        cme_time = self.ui.CME_time.time().toString("hh:mm AP")

        # Combine date and time into a single string
        time_utc = f"{cme_date} {cme_time}"  # Example: "2000-01-01 00:00"

        # Convert the combined string into a datetime object
        # This variable convert the date time object to comatible format for the code.
        # converting string to date time object
        time_utc = datetime.strptime(time_utc, "%d-%m-%Y %I:%M %p")
        T0 = time_utc.timestamp()
        ic(T0)

        # Get values from spin boxes
        cme_dt0 = self.ui.CME_dt0.value()
        ic(cme_dt0)
        ic()
        cme_dt0 = cme_dt0*60.0
        ic(cme_dt0)
        r0 = self.ui.R0.value()
        r0 = ((r0*u.R_sun).to(u.km)).value
        dr0 = self.ui.dR0.value()
        dr0 = ((dr0*u.R_sun).to(u.km)).value
        v0 = self.ui.V0.value()
        dv0 = self.ui.dV0.value()
        w_in = self.ui.w_in.value()
        dw = self.ui.dW.value()
        gamma = self.ui.Gamma.value()
        dgamma = self.ui.dGamma.value()

        # Get value from combo box and updating for calculation
        target = self.ui.Target.currentText()
        _, R_target = f.position(target, time_utc)

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
                self.reset_1D()  # Reset inputs if "OK" is clicked
            return  # Exit the function to prevent further calculations
        r1 = ((R_target*u.au).to(u.km)).value
        
       
        ic(r0)
        ic(dr0)
        ic(r1)

        # Doing DBM calculation for 1D input
        if pdbm_choice == "Yes" and auto_choice == "Yes":

            # Number of Ensembles
            N = 10000

            if wind_type == 'Slow':
                wind_array = np.clip(st.norm.rvs(
                    370.530847, 88.585045/3.0, size=N), 1, 1000)
                gamma_array = np.clip(st.lognorm.rvs(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08, size=N), 1.0e-08, 3.0e-7)

            else:
                wind_array = np.clip(st.norm.rvs(
                    579.057905, 67.870776/3.0, size=N), 1, 1000)
                gamma_array = np.clip(st.lognorm.rvs(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08, size=N), 1.0e-08, 3.0e-7)

            T_array, V_array, V0_array = f.PDBM(
                r0, dr0, r1, v0, dv0, gamma_array, wind_array, cme_dt0, N)

            T_mean = np.nanmean(T_array)
            T_std = np.nanstd(T_array)
            V_mean = np.nanmean(V_array)
            V_std = np.nanstd(V_array)
            W_median = np.nanmedian(wind_array)
            G_median = np.nanmedian(gamma_array)
            t_arrival = T0+(np.nanmedian(T_array)*3600)
            t_arrival_UTC = datetime.fromtimestamp(
                t_arrival).strftime("%Y-%m-%d %H:%M")
            tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
            T_arrival_UTC = datetime.strptime(t_arrival_UTC, "%Y-%m-%d %H:%M")

            output_text = (f"* CME arrival forecast for <b>{target}</b>.<br>"
                           f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                           f"* Mean transit time of CME is <b>{T_mean:.3f} hrs</b> and uncertainty is <b>{T_std:.3f} hrs</b>.<br>"
                           f"* Mean impact speed of CME is <b>{V_mean:.3f} km/s</b> and uncertainty is <b>{V_std:.3f} km/s</b>.<br>"
                           f"* Median value of drag parameter γ is <b>{G_median/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                           f"* Median Value of solar wind speed w is <b>{W_median:.2f} km/s</b>.</br> ")

            # Making plot
            plot = f.PDBM_RVT_plot(
                time_utc, T_array, r0, V0_array, gamma_array, wind_array, R_target, tdate)
            T_PDF_plot = f.TT_plot(T_array)
            V_PDF_Plot = f.V_plot(V_array)

            # Showing Output
            self.ui.DBM_Result_1D.show()
            self.ui.DBM_Result_1D.clear()
            self.ui.DBM_Result_1D.setText(output_text)
            self.ui.DBM_Result_1D.adjustSize()
            self.ui.DBM_Result_1D.setFixedWidth(551)

            # Showing Plot
            pixmap = QPixmap()
            pixmap.loadFromData(plot.read())
            self.ui.RVT_plot.show()
            self.ui.RVT_plot.clear()
            self.ui.RVT_plot.setGeometry(10, 440, 500, 381)
            self.ui.RVT_plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.RVT_plot.setScaledContents(True)

            pixmap = QPixmap()
            pixmap.loadFromData(T_PDF_plot.read())
            self.ui.TT_Plot.show()
            self.ui.TT_Plot.clear()
            self.ui.TT_Plot.setGeometry(515, 440, 400, 381)
            self.ui.TT_Plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.TT_Plot.setScaledContents(True)

            pixmap = QPixmap()
            pixmap.loadFromData(V_PDF_Plot.read())
            self.ui.V_plot.show()
            self.ui.V_plot.clear()
            self.ui.V_plot.setGeometry(920, 440, 400, 381)
            self.ui.V_plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.V_plot.setScaledContents(True)

        elif pdbm_choice == "Yes" and auto_choice == "No":
            # Number of Ensembles
            N = 10000

            wind_array = np.random.normal(w_in, dw/3.0, N)
            gamma_array = np.clip(np.random.normal(
                gamma*1.0e-7, dgamma*1.0e-7/3.0, N), 1.0e-15, 3.0e-2)

            T_array, V_array, V0_array = f.PDBM(
                r0, dr0, r1, v0, dv0, gamma_array, wind_array, cme_dt0, N)

            T_mean = np.nanmean(T_array)
            T_std = np.nanstd(T_array)
            V_mean = np.nanmean(V_array)
            V_std = np.nanstd(V_array)
            W_median = np.nanmedian(wind_array)
            G_median = np.nanmedian(gamma_array)
            t_arrival = T0+(np.nanmedian(T_array)*3600)
            t_arrival_UTC = datetime.fromtimestamp(
                t_arrival).strftime("%Y-%m-%d %H:%M")
            tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
            T_arrival_UTC = datetime.strptime(t_arrival_UTC, "%Y-%m-%d %H:%M")

            output_text = (f"* CME arrival forecast for <b>{target}</b>.<br>"
                           f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                           f"* Mean transit time of CME is <b>{T_mean:.3f} hrs</b> and uncertainty is <b>{T_std:.3f} hrs</b>.<br>"
                           f"* Mean impact speed of CME is <b>{V_mean:.3f} km/s</b> and uncertainty is <b>{V_std:.3f} km/s</b>.<br>"
                           f"* Median value of drag parameter γ is <b>{G_median/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                           f"* Median Value of solar wind speed w is <b>{W_median:.2f} km/s</b>.</br> ")

            # Making plot
            plot = f.PDBM_RVT_plot(
                time_utc, T_array, r0, V0_array, gamma_array, wind_array, R_target, tdate)
            T_PDF_plot = f.TT_plot(T_array)
            V_PDF_Plot = f.V_plot(V_array)

            # Showing Output
            self.ui.DBM_Result_1D.show()
            self.ui.DBM_Result_1D.clear()
            self.ui.DBM_Result_1D.setText(output_text)
            self.ui.DBM_Result_1D.adjustSize()
            self.ui.DBM_Result_1D.setFixedWidth(551)

            # Showing Plot
            pixmap = QPixmap()
            pixmap.loadFromData(plot.read())
            self.ui.RVT_plot.show()
            self.ui.RVT_plot.clear()
            self.ui.RVT_plot.setGeometry(10, 440, 500, 381)
            self.ui.RVT_plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.RVT_plot.setScaledContents(True)

            pixmap = QPixmap()
            pixmap.loadFromData(T_PDF_plot.read())
            self.ui.TT_Plot.show()
            self.ui.TT_Plot.clear()
            self.ui.TT_Plot.setGeometry(515, 440, 400, 381)
            self.ui.TT_Plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.TT_Plot.setScaledContents(True)

            pixmap = QPixmap()
            pixmap.loadFromData(V_PDF_Plot.read())
            self.ui.V_plot.show()
            self.ui.V_plot.clear()
            self.ui.V_plot.setGeometry(920, 440, 400, 381)
            self.ui.V_plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.V_plot.setScaledContents(True)

        elif pdbm_choice == "No" and auto_choice == "Yes":
            if wind_type == 'Slow':
                W_1D = st.norm.mean(
                    370.530847, 88.585045)
                gamma_1D = st.lognorm.median(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08)
            else:
                W_1D = st.norm.mean(
                    579.057905, 67.870776)
                gamma_1D = st.lognorm.median(
                    0.6518007540612114, -2.2727287735377082e-08, 9.425812152200486e-08)

            TT, VT = f.DBM(r0, r1, v0, gamma_1D, W_1D)
            t_arrival = T0+(TT*3600)
            t_arrival_UTC = datetime.fromtimestamp(
                t_arrival).strftime("%Y-%m-%d %H:%M")
            tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
            T_arrival_UTC = datetime.strptime(t_arrival_UTC, "%Y-%m-%d %H:%M")

            output_text = (f"* CME arrival forecast for <b>{target}</b>.<br>"
                           f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                           f"* Mean transit time of CME is <b>{TT:.3f} hrs</b>.<br>"
                           f"* Mean impact speed of CME is <b>{VT:.3f} km/s</b>.<br>"
                           f"* Median value of drag parameter γ is <b>{gamma_1D/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                           f"* Median Value of solar wind speed w is <b>{W_1D:.2f} km/s</b>.</br> ")

            # Making Plots
            plot = f.DBM_RVT_plot(time_utc, TT, r0,
                                  v0, gamma_1D, W_1D, R_target, tdate)

            # Showing Output
            self.ui.DBM_Result_1D.show()
            self.ui.DBM_Result_1D.clear()
            self.ui.DBM_Result_1D.setText(output_text)
            self.ui.DBM_Result_1D.adjustSize()
            self.ui.DBM_Result_1D.setFixedWidth(551)

            # Showing Plot
            self.ui.TT_Plot.hide()
            self.ui.V_plot.hide()
            pixmap = QPixmap()
            pixmap.loadFromData(plot.read())
            self.ui.RVT_plot.show()
            self.ui.RVT_plot.clear()
            self.ui.RVT_plot.setGeometry(120, 400, 700, 400)
            self.ui.RVT_plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.RVT_plot.setScaledContents(True)

        elif pdbm_choice == "No" and auto_choice == "No":
            W_1D = w_in
            gamma_1D = gamma*1.0e-7
            TT, VT = f.DBM(r0, r1, v0, gamma_1D, W_1D)
            t_arrival = T0+(TT*3600)
            t_arrival_UTC = datetime.fromtimestamp(
                t_arrival).strftime("%Y-%m-%d %H:%M")
            tdate = datetime.strptime(t_arrival_UTC, '%Y-%m-%d %H:%M')
            T_arrival_UTC = datetime.strptime(t_arrival_UTC, "%Y-%m-%d %H:%M")

            output_text = (f"* CME arrival forecast for <b>{target}</b>.<br>"
                           f"* CME will arrive at <b>{R_target:.3f} AU</b> on (date and time) <b>{T_arrival_UTC}</b>.<br>"
                           f"* Mean transit time of CME is <b>{TT:.3f} hrs</b>.<br>"
                           f"* Mean impact speed of CME is <b>{VT:.3f} km/s</b>.<br>"
                           f"* Median value of drag parameter γ is <b>{gamma_1D/1e-7:.3f}e-7 km<sup>-1</sup></b>.<br>"
                           f"* Median Value of solar wind speed w is <b>{W_1D:.2f} km/s</b>.</br> ")

            # Making Plots
            plot = f.DBM_RVT_plot(time_utc, TT, r0,
                                  v0, gamma_1D, W_1D, R_target, tdate)

            # Showing Output
            self.ui.DBM_Result_1D.show()
            self.ui.DBM_Result_1D.clear()
            self.ui.DBM_Result_1D.setText(output_text)
            self.ui.DBM_Result_1D.adjustSize()
            self.ui.DBM_Result_1D.setFixedWidth(551)

            # Showing Plot
            self.ui.TT_Plot.hide()
            self.ui.V_plot.hide()
            pixmap = QPixmap()
            pixmap.loadFromData(plot.read())
            self.ui.RVT_plot.show()
            self.ui.RVT_plot.clear()
            self.ui.RVT_plot.setGeometry(120, 400, 700, 400)
            self.ui.RVT_plot.setPixmap(pixmap)
            # Ensure the plot scales to fit the QLabel
            self.ui.RVT_plot.setScaledContents(True)

  