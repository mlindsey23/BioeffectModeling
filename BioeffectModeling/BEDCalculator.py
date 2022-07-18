'''
Convert Absorbed Dose to Biological Effective Dose (BED)

Input -> RTDOSE, RTSTRUCT, CT
Output -> RTDOSE
'''

import numpy as np
from datetime import datetime
from DICOM_RT import DicomPatient as dcmpat
from MIRDCalculation_BED.BioeffectModeling.ROI_Values import *


# Note: Only includes branches for Normal Liver, Normal Lungs, and other
# That's all I could find in the RTSTRUCT files I tested, but simple to add more

def GetBEDinDICOM(basepath, dosepath, structpath, HighestVoxelValue, ct_path = '', nm_path = ''):
    if ct_path == '':
        ctpath = basepath + '/CT/'
    else:
        ctpath = ct_path
    if nm_path == '':
        nmpath = basepath + '/NM/'
    else:
        nmpath = nm_path
    
    calc = BioeffectCalculator(basepath, ctpath, nmpath, structpath, dosepath)
    
    if calc.ctObject.img3D.shape[2] < calc.activityObject.img3D.shape[2]:
        BEDimg3D = np.zeros((calc.activityObject.img3D.shape[0],calc.activityObject.img3D.shape[1],calc.ctObject.img3D.shape[2]))
    else :
        BEDimg3D = np.zeros(calc.activityObject.img3D.shape)

    calc.BEDCalculator(BEDimg3D)
    BEDimg3D = BEDimg3D / np.max(BEDimg3D) * HighestVoxelValue
    calc.WriteRTDoseCT(BEDimg3D)
    print("End BED Calculation.")

    
class BioeffectCalculator(dcmpat.PatientCT, dcmpat.Patient3DActivity):
    def __init__(self, basepath, ctpath, nmpath, structpath, dosepath):
        self.patientObject = dcmpat.DicomPatient(basepath)
        self.ctObject = dcmpat.PatientCT(ctpath)
        self.activityObject = dcmpat.Patient3DActivity(nmpath)
        self.activityObject.LoadRTDose(dosepath)
        self.ctObject.LoadStructures(structpath)

    def BEDCalculator(self, BEDimg3D):
            di = round(self.ctObject.img3D.shape[0]/self.activityObject.img3D.shape[0])
            dj = round(self.ctObject.img3D.shape[1]/self.activityObject.img3D.shape[1])
            for i in range(self.activityObject.img3D.shape[0]):
                prog = i/self.activityObject.img3D.shape[0]*100
                print("Calculating BED... (" + str(round(prog,1))+"%)")
                for j in range(self.activityObject.img3D.shape[1]): 
                    if self.ctObject.img3D.shape[2] < self.activityObject.img3D.shape[2]:
                        for k in range(self.ctObject.img3D.shape[2]):
                            if self.ctObject.structures3D['Liver'][i*di,j*dj,k] == True :
                                if self.ctObject.structures3D['All Tumors (Left Lobe)'][i*di,j*dj,k] == True or self.ctObject.structures3D['All Tumors (Right Lobe)'][i*di,j*dj,k] == True:
                                    Trep = Trep_Tumor
                                    AlphaBeta = AlphaBeta_TLiver
                                else: 
                                    Trep = Trep_Normal
                                    AlphaBeta = AlphaBeta_NLiver
                            elif self.ctObject.structures3D['Lung_L'][i*di,j*dj,k] == True or self.ctObject.structures3D['Lung_R'][i*di,j*dj,k] == True :
                                Trep = Trep_Normal
                                AlphaBeta = AlphaBeta_NLung
                            else :
                                Trep = Trep_Normal
                                AlphaBeta = AlphaBeta_Standard      
                            BEDimg3D[i,j,k] = self.activityObject.img3D[i,j,k] * (1 + (( self.activityObject.img3D[i,j,k] * Trep) / (AlphaBeta * (Trep + RadionuclideHalfLife))))
                    else:
                        for k in range(self.activityObject.img3D.shape[2]):
                            if self.ctObject.structures3D['Liver'][i*di,j*dj,k] == True :
                                if self.ctObject.structures3D['All Tumors (Left Lobe)'][i*di,j*dj,k] == True or self.ctObject.structures3D['All Tumors (Right Lobe)'][i*di,j*dj,k] == True:
                                    Trep = Trep_Tumor
                                    AlphaBeta = AlphaBeta_TLiver
                                else: 
                                    Trep = Trep_Normal
                                    AlphaBeta = AlphaBeta_NLiver
                            elif self.ctObject.structures3D['Lung_L'][i*di,j*dj,k] == True or self.ctObject.structures3D['Lung_R'][i*di,j*dj,k] == True :
                                Trep = Trep_Normal
                                AlphaBeta = AlphaBeta_NLung
                            else :
                                Trep = Trep_Normal
                                AlphaBeta = AlphaBeta_Standard
                            BEDimg3D[i,j,k] = self.activityObject.img3D[i,j,k] * (1 + (( self.activityObject.img3D[i,j,k] * Trep) / (AlphaBeta * (Trep + RadionuclideHalfLife))))
    
    def WriteRTDoseCT(self, BEDimg3D):
        self.ctObject.WriteRTDose(BEDimg3D, name = 'BEDCalculation_' + datetime.now().strftime("%m%d%y_%H%M%S") + '.dcm' ) 
