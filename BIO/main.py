from BioeffectModeling.BIO.EUBEDCalculator import *


### USER PARAMETERS ###

# 1. Path to DICOM files (str)
basepath = '/Users/mjlindsey/Documents/LiverPatients/Patient1/POSTTX/'

# 2. RTDOSE filename (str)
dosefile = 'TOPASDose.dcm'

# 3. Units for BED RTDOSE file (str) (Choose from "Gy/GBq", "Gy/mCi", or "Gy")
unit = "Gy"

# 4. ROIs for EUBED Calculation (str, list)
ROIList = ['Liver', 'Lung_L', 'Lung_R']

# 5. Print EUBED and EUD into txt file (bin)
createFile = True

# 6. Scale input dose file to a max value (float) (Note: Useful for low-dose TOPAS files) (optional)
maxVoxel = 50

### Main Script ###
calc = EUBEDCalculator(basepath, dosefile, unit, maxVoxel)
calc.BEDCalculator()    
calc.WriteRTDoseBED()
calc.EUBED(ROIList, createFile)
calc.EUD(ROIList, createFile)
