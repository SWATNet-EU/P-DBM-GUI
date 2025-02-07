# P-DBM
Probabilistic Drag Based Model (P-DBM) is probabilistic approach on Drag Based Model (DBM) to predict arrival time and speed of coronal mass ejection (CME) at any target in the heliosphere.
Detailed description of the models can be found in [Vršnak et al., 2013 (DBM)](https://doi.org/10.1007/s11207-012-0035-4) and [Napoletano et al., 2018 (P-DBM)](https://doi.org/10.1051/swsc/2018003).

In the following, how DBM and P-DBM calaculations can be performed using this graphical user interface (GUI). In the probabilitic approach, a probability distribution function is buit around the defined input parameters of the CME and DBM.  

## Setup
This GUI has been created with pyQT6 framework which is comaptible with a python-3.9.
The most recent LTS version of Ubuntu is equipped with a python-3.12; to avoid breaking an OS it recommanded to have working installation of Miniconda/Anaconda.

The following commands will create a virtual enviroment to run this gui on your local system.
```
conda create --name gui python=3.9.21
```
To activate the enviroment:
```
conda activate gui
```
Next, clone the repository with following command:
```
git clone https://github.com/astronish16/P-DBM-GUI.git
```
Move to the project's directory
```
cd P-DBM-GUI
```
Now, install all the necessary python packages to the virtual enviroment using the provided command:
```
pip install -r requirements.txt
```

After this GUI can be run on local machine by hitting
```
python3 gui_v1.1.py
```

## P-DBM-GUI
As an output of above command, The following window will be pop up to asking an input for the (P)-DBM calculations.

![1D-DBM](https://github.com/user-attachments/assets/ac88ec5d-fae3-43f4-92ff-872630c1b6f4)

![1D-PDBM-Output](https://github.com/user-attachments/assets/44e7f3fe-341d-4dbf-bc59-cc3cc52b1733)
