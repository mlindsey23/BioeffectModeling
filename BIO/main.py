from BioeffectModeling.BIO.EUBEDCalculator import *


### USER PARAMETERS ###

# 1. Path to DICOM files (str)
basepath = '/Users/mjlindsey/Documents/LiverPatients/Patient1/POSTTX/'

# 2. RTDOSE filename (str)
dosefile = 'MIRDDose.dcm'

# 3. Units for BED RTDOSE file (str) (Choose from "Gy/Gbq", "Gy/mCi", or "Gy")
unit = "Gy/GBq"

# 4. ROIs for EUBED Calculation (str, list)
ROIList = ['Liver', 'Lung_L', 'Lung_R']




### Main Script ###
calc = EUBEDCalculator(basepath, dosefile, unit)
calc.BEDCalculator()    
calc.WriteRTDoseBED()
calc.EUBED(ROIList, CreateFile)
