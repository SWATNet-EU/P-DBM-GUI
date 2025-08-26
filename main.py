import streamlit as st
import tab1
import tab2

# Set the page configuration
st.set_page_config(
    page_title="CME Arrival Time Calculator", layout="wide", page_icon="☀️"
)

# Main title
st.title("☀️ P-DBM GUI")

# Create tabs
dbm1, dbm2, us = st.tabs(["1D-DBM", "2D - DBM", "About Us"])

# ----- TAB 1 -----
with dbm1:
    tab1.display()

# ----- TAB 2 -----
with dbm2:
    tab2.display()

# ----- TAB 3 -----
with us:
    col1, col2, col3 = st.columns(3)

    # with col1:
    #     st.image("About_Us_logo/UNITOV.jpg", width=120)

    with col2:
        col201, col202, col203 = st.columns([1, 2, 1])
        with col201:
            st.image("About_Us_logo/UNITOV.jpg", width=150)

        with col202:
            st.image("About_Us_logo/SWATNet.png", width=300)

        with col203:
            st.image("About_Us_logo/UOS.png", width=150)

    # with col3:
    #     st.image("About_Us_logo/UOS.png", width=120)

    st.markdown(
        """
    <div style="text-align: center;">
        <h2 style="color:#11262f;">🖥️ <b>P-DBM GUI</b> – A User-Friendly Web App for CME Forecasting</h2>
    </div>

    <p style="font-size:16pt; text-align: justify;">
        Welcome to <strong>P-DBM GUI</strong>, an intuitive and interactive web application designed to estimate Coronal Mass Ejection (CME) arrival times using the <strong>Probabilistic Drag Based Model (P-DBM)</strong>.
        This tool aims to perform P-DBM simulation result and provide visulization to the space weather community.
        Developed under the <a href="https://swatnet.eu/project9/" target="_blank" style="text-decoration: none;"><strong>SWATNet Project-9</strong></a>, this app is part of our commitment to make space weather tools more accessible and scientifically robust.
    </p>

    <div style="background-color:#d3fbd8; padding: 10px; border-left: 10px solid #1f77b4; margin-top: 20pt;font-size:16pt">
        <strong>📌 Acknowledgement:</strong> <br>
        SWATNet project is funded by the European Union's Horizon 2020 programme under the Marie Skłodowska-Curie Grant Agreement No. <strong>955620</strong>.
    </div>

    <div style="background-color:#f9f9f9; padding: 10px; border-left: 10px solid #1f77b4; margin-top: 20pt">
        <strong>📑 Citation:</strong> Please cite the following articles if you use P-DBM GUI in your publication.<br>
        (1) A catalogue of observed geo-effective CME/ICME characteristics.
            Ronish Mugatwala, Simone Chierichini, Gregoire Francisco, Gianluca Napoletano, Raffaello Foldes, Luca Giovannelli, Giancarlo De Gasperis, Enrico Camporeale, Robertus Erdélyi and Dario Del Moro.
            J. Space Weather Space Clim., 14 (2024) 6.
            DOI: <a href="https://doi.org/10.1051/swsc/2024004">https://doi.org/10.1051/swsc/2024004</a>
        
    </div>

    <p style="text-align: center; margin-top: 30pt;">
        🙏 <em>Thank you for using our tool!</em>
    </p>
    """,
        unsafe_allow_html=True,
    )
