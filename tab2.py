import streamlit as st
import matplotlib.pyplot as plt
from datetime import datetime
from icecream import ic
from calculation_1D import DBM1D
from calculation_2D import DBM2D
import dbm_functions as dbm
from PIL import Image
from io import BytesIO

# By default st.date_input doesn't allow to provide date before 2013.
# so we create a variable to store minimum date for our use case.
min_date = datetime(1995, 1, 1)  # Example: January 1, 2000


def display():
    col1, col2, col3 = st.columns([0.25, 0.2, 0.4], border=True)

    # ? The following values are incorporated to resolve unbound error.
    # ? This values are not used in calculation at all.
    dt0, dr0, dv0, dw, dgamma, w, gamma, solar_wind, domega, dphi = (
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        "Slow",
        None,
        None,
    )  # This will be ignored

    with col1:
        pdbm_col1, pdbm_col2 = st.columns([0.3, 0.3])
        with pdbm_col1:
            st.markdown("P-DBM Calculation:", unsafe_allow_html=True)
        with pdbm_col2:
            pdbm = st.radio(
                "",
                ("Yes", "No"),
                horizontal=True,
                label_visibility="collapsed",
                key="pdbm_2d",
            )

        st.markdown("##### <u> CME Input Parameters </u>", unsafe_allow_html=True)

        with st.container():
            # Row 1: CME Date
            row1_col1, row1_col2, row1_col3 = st.columns([1.2, 1, 1])
            with row1_col1:
                st.markdown("CME Date at $R_0$:")
            with row1_col2:
                cme_date = st.date_input(
                    "",
                    key="cme_date_2d",
                    label_visibility="collapsed",
                    value=datetime.today().date(),
                    help="Select date for the CME Observation using pop-up calender. You can also write it in YYYY/MM/DD format.",
                    min_value=min_date,
                )

            # Row 2: CME Time
            row2_col1, row2_col2, row2_col3, row2_col4 = st.columns(
                [0.85, 0.45, 0.45, 0.45]
            )
            with row2_col1:
                st.markdown("CME Time at $R_0$:")
            with row2_col2:
                cme_hour = st.selectbox(
                    "",
                    options=list(range(0, 24)),
                    index=12,
                    help="CME hour",
                    key="cme_hour_2d",
                )
            with row2_col3:
                cme_min = st.selectbox(
                    "",
                    options=list(range(0, 60)),
                    index=0,
                    help="CME minute",
                    key="cme_ninute_2d",
                )

            if pdbm == "Yes":
                with row2_col4:
                    dt0 = st.selectbox(
                        "",
                        options=list(range(0, 120)),
                        index=60,
                        help="Uncertanity in CME time (in minutes)",
                        key="cme_dt0_2d",
                    )

            # * CME date used for calculation
            cme_date_time = datetime.combine(cme_date, datetime.min.time()).replace(
                hour=cme_hour, minute=cme_min
            )

            # Row 3: CME Initial Position
            row3_col1, row3_col2, row3_col3 = st.columns([1.2, 1, 1])
            with row3_col1:
                st.markdown("CME Initial Position $R_0(R_\odot)$:")
            with row3_col2:
                r0 = st.number_input(
                    "",
                    key="r0_2d",
                    min_value=20.0,
                    value=20.0,
                    label_visibility="collapsed",
                )
            if pdbm == "Yes":
                with row3_col3:
                    dr0 = st.number_input(
                        "",
                        min_value=0.1,
                        value=1.0,
                        label_visibility="collapsed",
                        help="Uncertanity in R0 (R_sun)",
                        key="dr0_2d",
                    )

            # Row 4: CME Speed
            row4_col1, row4_col2, row4_col3 = st.columns([1.2, 1, 1])
            with row4_col1:
                st.markdown("Speed of CME at $R_0$ - $V_0$ (km/s):")
            with row4_col2:
                v0 = st.number_input(
                    "", key="v0_2d", value=1000.0, label_visibility="collapsed"
                )
            if pdbm == "Yes":
                with row4_col3:
                    dv0 = st.number_input(
                        "",
                        min_value=0.0,
                        value=100.0,
                        label_visibility="collapsed",
                        key="dv0_2d",
                    )

            # Row 5: Target
            row5_col1, row5_col2 = st.columns([1, 2])
            with row5_col1:
                st.markdown("Target:")
            with row5_col2:
                target_2d = st.selectbox(
                    "",
                    [
                        "Mercury",
                        "Venus",
                        "Earth",
                        "Mars",
                        "Jupiter",
                        "Saturn",
                        "Messenger",
                        "VEX",
                        "PSP",
                        "SolO",
                        "BepiCol",
                        "Spitzer",
                        "Wind",
                        "ST-A",
                        "ST-B",
                        "Kepler",
                        "MSL",
                        "Maven",
                        "Juno",
                    ],
                    key="target_2d",
                    label_visibility="collapsed",
                )

    with col2:
        st.markdown("##### <u>CME Morphology Input</u>", unsafe_allow_html=True)
        with st.container():
            row6_col1, row6_col2 = st.columns([0.3, 0.7])
            # Row-6 Cone type
            with row6_col1:
                st.markdown("Cone Type:")
            with row6_col2:
                cone_type = st.radio(
                    "",
                    ("Ice-Cream Cone", "Tangential Cone", "Concentric Cone"),
                    horizontal=True,
                    key="cone_type_2d",
                    label_visibility="collapsed",
                )

            # Row-7 Kinematic approach
            row7_col1, row7_col2 = st.columns([0.3, 0.7])
            with row7_col1:
                st.markdown("Kinematic Type:")
            with row7_col2:
                kinematic_type = st.radio(
                    "",
                    ("Self-Similar Expansion", "Flattening Cone Evolution"),
                    horizontal=True,
                    key="kinematic_type_2d",
                    label_visibility="collapsed",
                )

            # Row-8 Half-Angular Width
            row8_col1, row8_col2, row8_col3 = st.columns([1.2, 1, 1])
            with row8_col1:
                st.markdown("Half Angular Width $\omega$:")
            with row8_col2:
                omega = st.number_input(
                    "",
                    key="omega_2d",
                    value=60,
                    max_value=90,
                    label_visibility="collapsed",
                )

            if pdbm == "Yes":
                with row8_col3:
                    domega = st.number_input(
                        "",
                        key="domega_2d",
                        value=10,
                        label_visibility="collapsed",
                    )

            # Row-9 CME Propagation direction
            row9_col1, row9_col2, row9_col3 = st.columns([1.2, 1, 1])
            with row9_col1:
                st.markdown("CME Propagation Direction $\Phi_{CME}$:")
            with row9_col2:
                phi_cme = st.number_input(
                    "",
                    key="phi_2d",
                    value=0,
                    max_value=360,
                    label_visibility="collapsed",
                )

            if pdbm == "Yes":
                with row9_col3:
                    dphi = st.number_input(
                        "",
                        key="dphi_2d",
                        value=10,
                        label_visibility="collapsed",
                    )
            st.markdown("---")
            st.markdown("##### <u>DBM Input Parameters</u>", unsafe_allow_html=True)

            # Row-1 Auto Choice for DBM
            row1_col1, row1_col2 = st.columns([0.5, 0.5])
            with row1_col1:
                st.markdown("Auto Choice:")
            with row1_col2:
                auto_dbm = st.radio(
                    "",
                    ("Yes", "No"),
                    horizontal=True,
                    key="autodbm_2d",
                    label_visibility="collapsed",
                )

            # Row-2 Solar wind type choice
            if auto_dbm == "Yes":
                row2_col1, row2_col2 = st.columns([0.5, 0.5])
                with row2_col1:
                    st.markdown("Solar Wind Type:")
                with row2_col2:
                    solar_wind = st.radio(
                        "",
                        ("Slow", "Fast"),
                        horizontal=True,
                        key="solarwind_2d",
                        label_visibility="collapsed",
                    )

            # Sola wind speed manual
            if auto_dbm == "No":
                row3_col1, row3_col2, row3_col3 = st.columns([0.5, 0.25, 0.25])

                with row3_col1:
                    st.markdown("Ambient Solar Wind Speed $w$ (km/s):")
                with row3_col2:
                    w = st.number_input(
                        "",
                        value=400.0,
                        min_value=0.0,
                        label_visibility="collapsed",
                        key="2D_w",
                    )

                if pdbm == "Yes":
                    with row3_col3:
                        dw = st.number_input(
                            "",
                            value=30.0,
                            min_value=0.0,
                            label_visibility="collapsed",
                            help="Uncertanity in solar wind speed",
                            key="2D_dw",
                        )
                # drag parameter manual
                row4_col1, row4_col2, row4_col3 = st.columns([0.5, 0.25, 0.25])

                with row4_col1:
                    st.markdown("Drag Parameter $\Gamma$ (x 10$^{-7}$ km$^{-1}$):")
                with row4_col2:
                    gamma = st.number_input(
                        "",
                        value=0.2,
                        min_value=0.001,
                        label_visibility="collapsed",
                        key="2D_gamma",
                    )

                if pdbm == "Yes":
                    with row4_col3:
                        dgamma = st.number_input(
                            "",
                            value=0.1,
                            min_value=0.0001,
                            label_visibility="collapsed",
                            help="Uncertanity in drag parameter",
                            key="2D_dgamma",
                        )

            st.markdown("---")
            calculate_2d = st.button("Calculate", key="2d_calculation")

    if calculate_2d:
        with st.spinner("Performing P-DBM Calculation"):
            if pdbm == "Yes":
                pdbm = True
            else:
                pdbm = False

            if auto_dbm == "Yes":
                auto_dbm = True
            else:
                auto_dbm = False

            dbm2d = DBM2D(
                time_utc=cme_date_time,
                r0=r0,
                v0=v0,
                Omega=omega,
                Phi_CME=phi_cme,
                cone_type=cone_type,
                Kinematic=kinematic_type,
                target_name=target_2d,
                P_DBM=pdbm,
                auto_dbm=auto_dbm,
                wind_type=solar_wind,
                w=w,
                gamma=gamma,
                dt=dt0,
                dr0=dr0,
                dv0=dv0,
                dw=dw,
                dgamma=dgamma,
                domega=domega,
                dphi_cme=dphi,
            )
            predictions = dbm2d.P_DBM_run()
            ic(predictions)

            with col3:
                if pdbm:
                    st.markdown(" #### (P)-DBM Results: 📊", unsafe_allow_html=True)
                    st.markdown(
                        f"""
                    **CME arrival date and time:** {predictions["Arrival_time"]} with **Probability of Arrival:** {predictions["Probability of Arrival"]:.3f}   
                    **Transit time of CME:** {predictions["Transit_time_mean"]:.2f} $\pm$ {predictions["Transit_time_std"]:.2f} hr.  
                    **Transit time spread (90%):** {predictions["T_10"]:.2f} - {predictions["T_90"]:.2f} hr.  
                    **Arrival speed of CME at {predictions["Target"]}:** {predictions["Arrival_speed_mean"]:.2f} $\pm$ {predictions["Arrival_speed_std"]:.2f} km/s.  
                    **Arrival speed spread (90%):** {predictions["V_10"]:.2f} - {predictions["V_90"]:.2f} km/s.  
                    **CME travel distance:** {predictions["Travel_distance"]:.3f} AU.  
                    CME launched from {r0:.2f} $R_\odot$ with **initial speed** {v0:.2f} km/s on {predictions["Initial_time"]}."""
                    )

                    st.markdown("---")
                    st.markdown("#### (P)-DBM Model Parameterrs:")
                    st.markdown(
                        f"""
                    **Ambient Solar wind speed w** = {predictions["w_median"]:.2f} km/s.  
                    **Drag parameter $\gamma$** : {predictions["gamma_median"]:2E} km$^{-1}$. """
                    )
                else:
                    st.markdown(" #### (P)-DBM Results:", unsafe_allow_html=True)
                    st.markdown(
                        f"""
                    **CME arrival date and time:** {predictions["Arrival_time"]}.  
                    **Transit time of CME:** {predictions["Transit_time_mean"]:.2f} hr.  
                    **Arrival speed of CME at {predictions["Target"]}:** {predictions["Arrival_speed_mean"]:.2f} km/s.  
                    **CME travel distance:** {predictions["Travel_distance"]:.3f} AU.  
                    CME launched from {r0:.2f} $R_\odot$ with **initial speed** {v0:.2f} km/s on {predictions["Initial_time"]}."""
                    )

                    st.markdown("---")
                    st.markdown("#### (P)-DBM Model Parameterrs:")
                    st.markdown(
                        f"""
                    **Ambient Solar wind speed w** = {predictions["w_median"]:.2f} km/s.  
                    **Drag parameter $\gamma$** : {predictions["gamma_median"]:2E} km$^{-1}$. """
                    )
                CME = predictions["Heliosphere"]
                CME.seek(0)
                img = Image.open(CME)
                st.image(img)

        row2_c1, row2_c2, row2_c3 = st.columns([0.6, 0.4, 0.4], border=True)
        with row2_c1:
            RVT = predictions["RVT_plot"]
            RVT.seek(0)
            img = Image.open(RVT)
            st.image(img)

        if pdbm:
            with row2_c2:
                TT = predictions["TT_distribution"]
                TT.seek(0)
                img = Image.open(TT)
                st.image(img)
            with row2_c3:
                VV = predictions["V1_distribution"]
                VV.seek(0)
                img = Image.open(VV)
                st.image(img)
