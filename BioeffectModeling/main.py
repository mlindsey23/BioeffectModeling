from MIRDCalculation_BED.BioeffectModeling.BEDCalculator import *

### USER PARAMETERS ###

# 1. Path to DICOM files
basepath = '/Users/mjlindsey/Documents/Liver Patients/Patient12/'

# 2. Path to RTDOSE file
dosepath = '/Users/mjlindsey/Documents/Liver Patients/Patient12/DoseOnCTGrid.dcm'

# 3. Path to RTSTRUCT file
structpath = '/Users/mjlindsey/Documents/Liver Patients/Patient12/RTSTRUCT/2.16.840.1.114362.1.11972228.22981817573.591918584.369.1694.dcm'

# 4. Highest Voxel Value, to scale the output (NEEDS UPDATE)
HighestVoxelValue = 60

# 5. (optional) Alternate path to CT files (Default: ~/basepath/CT/ )
ctpath = ''

# 6. (optional) Alternate path to NM files (Default: ~/basepath/NM/ )
nmpath = ''



### MAIN ###

GetBEDinDICOM(basepath, dosepath, structpath, HighestVoxelValue, ctpath, nmpath)
