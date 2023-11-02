# Saq-aquifer-hydrological-model
CWatM model bespoke for the Saq aquifer. 

For standard CWatM: https://github.com/iiasa/CWatM.

Community Water Model (CWatM) is an open-source hydrological model simulating the water cycle daily at global and local levels, historically and into the future, maintained by IIASA's Water Security group. CWatM assesses water supply, demand, and environmental needs, including water management and human influence within the water cycle. CWatM includes an accounting of how future water demands will evolve in response to socioeconomic change and how water availability will change in response to climate and management.

Run with command 
python Saq-aquifer-hydrological-model\run_cwatm.py settings\settings-file.ini -l

Settings files relate to specific experiments. To tailor these for one's own system, the following are to be updated, including output folders and climate input:

- References to the Modflow folder (not Saudi_Arabia/Modflow) should point to Saq-aquifer-hydrological-model/Modflow. This comes compressed.
- References to the Saudi_Arabia folder should point to Saq-aquifer-hydrological-model/Saudi_Arabia.
- References to cwatm_input_5min2 refer to the cwatm_input_5min.zip folder on the IIASA Water FTP server: ftp://rcwatm:Water1090@ftp.iiasa.ac.at/
- metaNetcdfFile = Saq-aquifer-hydrological-model/metaNetcdf.xml


Note, Saq-aquifer-hydrological-model/Saudi_Arabia/Modflow is different than Saq-aquifer-hydrological-model/Modflow.
