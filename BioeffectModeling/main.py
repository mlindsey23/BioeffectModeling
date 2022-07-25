from MIRDCalculation_BED.BioeffectModeling.EUBEDCalculator import *


### USER PARAMETERS ###

# 1. Path to DICOM files (str)
basepath = '/Users/mjlindsey/Documents/LiverPatients/Patient12/'

# 2. RTDOSE filename (str)
dosefile = 'MIRDDose.dcm'

# 3. ROIs for EUBED Calculation (str, list)
ROIList = ['Liver', 'Lung_L', 'Lung_R']

# 4. Bins for DVH calculation (for faster processing) (posint)
bins = 2000

### Biological Effective Dose ###
calc = BioeffectCalculator(basepath, dosefile)
calc.BEDCalculator()    
calc.WriteRTDoseBED()
print("End BED Calculation.")

### EUBED ###
calc.EUBEDCalculator(ROIList)

### DVH Curves ###
dvhcalc = DVH(basepath, dosefile)
dvhcalc.DVHCalculator(ROIList, bins)
dvhcalc.PlotDVHCurves(ROIList, bins)
