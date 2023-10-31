# Saq-aquifer-hydrological-model
CWatM model bespoke for the Saq aquifer

Community Water Model (CWatM) is an open-source hydrological model simulating the water cycle daily at global and local levels, historically and into the future, maintained by IIASA's Water Security group. CWatM assesses water supply, demand, and environmental needs, including water management and human influence within the water cycle. CWatM includes an accounting of how future water demands will evolve in response to socioeconomic change and how water availability will change in response to climate and management.

Run with command 
python Saq-aquifer-hydrological-model\run_cwatm.py settings\settings-file.ini -l

To be updated in the settings files:

metaNetcdfFile = Saq-aquifer-hydrological-model/metaNetcdf.xml

References to the Modflow folder (not Saudi_Arabia/Modflow), now contained inside Saq-aquifer-hydrological-model
path_mf6dll = P:\watmodel\Modflow\mf6

References to the Saudi_Arabia folder, now contained inside Saq-aquifer-hydrological-model, for example
MaskMap = P:/watmodel/CWATM/Saudi_Arabia/mask_saq.tif

Saq-aquifer-hydrological-model/Saudi_Arabia/Modflow (referred to here) is different than Saq-aquifer-hydrological-model/Modflow (referred to in the aforementioned point)