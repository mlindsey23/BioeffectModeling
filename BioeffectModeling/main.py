from MIRDCalculation_BED.BioeffectModeling.EUBEDCalculator import *


### USER PARAMETERS ###

# 1. Path to DICOM files (str)
basepath = '/Users/mjlindsey/Documents/LiverPatients/Patient1/POSTTX/'

# 2. RTDOSE filename (str)
dosefile = 'MIRDDose.dcm'

# 3. ROIs for EUBED Calculation (str, list)
ROIList = ['Liver', 'Lung_L', 'Lung_R']

# 4. Create txt file of EUBED values (bin)
CreateFile = True

# 4. Bins for DVH calculation, to tune precision and processing time (posint)
bins = 2000

### Biological Effective Dose ###
calc = BioeffectCalculator(basepath, dosefile)
calc.BEDCalculator()    
calc.WriteRTDoseBED()
print("End BED Calculation.")

### EUBED ###
calc.EUBEDCalculator(ROIList, CreateFile)

### DVH Curves ###
dvhcalc = DVH(basepath, dosefile)
dvhcalc.DVHCalculator(ROIList, bins)
dvhcalc.PlotDVHCurves(ROIList, bins)
