'''
Convert Absorbed Dose to Biological Effective Dose (BED)
Input -> RTDOSE, RTSTRUCT, CT
Output -> RTDOSE
'''

import math
import os
import pydicom
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from BioeffectModeling.DICOM_RT import DicomPatient as dcmpat
from BioeffectModeling.BIO.ROI_Values import *
    
class EUBEDCalculator(dcmpat.PatientCT):
    def __init__(self, basepath, dosefile, unit = "Gy/GBq"):
        self.unit = unit
        ctpath = basepath + '/CT/'
        dosepath = basepath + dosefile
        dosefile_full = os.path.basename(dosefile)
        dosefile_split = dosefile_full.split('.')[0]
        self.basepath = basepath
        self.dosefilename = dosefile_split
        self.patientObject = dcmpat.DicomPatient(basepath)
        self.patientObject.dcmFileChosen = pydicom.dcmread(dosepath)
        self.ctObject = dcmpat.PatientCT(ctpath)
        self.ctObject.LoadRTDose(dosepath)
        try:
            structfile = os.listdir(basepath + '/RTSTRUCT/')
            structpath = basepath + '/RTSTRUCT/' + structfile[0]
            self.ctObject.LoadStructures(structpath)
        except Exception as e: 
            structfile = os.listdir(basepath + '/RTSTRUCT_LUNGSANDLIVER/')
            structpath = basepath + '/RTSTRUCT_LUNGSANDLIVER/' + structfile[0]
            self.ctObject.LoadStructures(structpath)
            print('ERROR: Could not load complete RTSTRUCT. CODE:', e)
            print('RTSTRUCT_LUNGSANDLIVER loaded instead.')
        self.ROIs = list(self.ctObject.structures3D.keys())
        print("ROI's identified:", self.ROIs)
        self.TUMORS = []
        for STRUCT in self.ROIs:
            if 'Tumor' in STRUCT:
                self.TUMORS.append(STRUCT)
        print("Tumor STRUCTS identified:", self.TUMORS)
        self.BEDimg3D = np.zeros(self.ctObject.img3D.shape)
        self.ConvertDoseUnits()

    def BEDCalculator(self):
        for i in range(self.ctObject.quantitiesOfInterest[0].array.shape[0]):
            if (i % 20) == 0:
                prog = i/self.ctObject.quantitiesOfInterest[0].array.shape[0]*100
                print("Calculating BED... (" + str(round(prog,1))+"%)")
            else:
                pass
            for j in range(self.ctObject.quantitiesOfInterest[0].array.shape[1]):
                for k in range(self.ctObject.quantitiesOfInterest[0].array.shape[2]):
                    if self.ctObject.structures3D['Liver'][i,j,k] == True :
                        for t in self.TUMORS:
                            if self.ctObject.structures3D[t][i,j,k] == True:
                                Trep = Trep_Tumor
                                AlphaBeta = AlphaBeta_TLiver
                                break
                            else:
                                Trep = Trep_Normal
                                AlphaBeta = AlphaBeta_NLiver
                    elif self.ctObject.structures3D['Lung_L'][i,j,k] == True or self.ctObject.structures3D['Lung_R'][i,j,k] == True :
                        Trep = Trep_Normal
                        AlphaBeta = AlphaBeta_NLung
                    else :
                        Trep = Trep_Normal
                        AlphaBeta = AlphaBeta_Standard
                    self.BEDimg3D[i,j,k] = self.ctObject.quantitiesOfInterest[0].array[i,j,k] * (1 + (( self.ctObject.quantitiesOfInterest[0].array[i,j,k] * Trep) / (AlphaBeta * (Trep + RadionuclideHalfLife))))
                                  
    def WriteRTDoseBED(self, seriesdescription = None):
        if seriesdescription == None:
            seriesdescription = 'BED_' + self.dosefilename
        name = 'BED_' + self.dosefilename + '.dcm'
        self.ctObject.WriteRTDose(self.BEDimg3D, self.basepath + name, self.unit, seriesdescription)

    def EUBED(self, ROIList, CreateFile):
        for r in ROIList:
            darr = np.zeros(self.BEDimg3D.shape)
            for i in range(self.BEDimg3D.shape[0]):
                if (i % 30) == 0:
                    prog = i/self.BEDimg3D.shape[0]*100
                    print("Calculating EUBED for " + r +  "... (" + str(round(prog,1))+"%)")
                for j in range(self.BEDimg3D.shape[1]):
                    for k in range(self.BEDimg3D.shape[2]):
                        if self.ctObject.structures3D[r][i,j,k] == True:
                            darr[i,j,k] = self.BEDimg3D[i,j,k]
            sum = 0
            if r not in self.TUMORS:
                Alpha_ROI = Alpha_Normal
            elif r in self.TUMORS:
                Alpha_ROI = Alpha_TLiver
            else:
                print('ROI type:' + r + ' not recognized')
                Alpha_ROI = Alpha_Normal
            for i in range(darr.shape[0]):
                for j in range(darr.shape[1]):
                    for k in range(darr.shape[2]):
                        if darr[i,j,k] > 0:
                            sum += math.exp(-Alpha_ROI * darr[i,j,k])
            N = (self.ctObject.structures3D[r]).sum()
            EUBED = -1/Alpha_Normal * np.log(sum / N)
            MEAN = np.sum(darr) / N
            RATIO = EUBED/MEAN
            if CreateFile == True:
                if ROIList.index(r) == 0:
                    f = open(self.basepath + 'EUBEDData_' + self.dosefilename + '.txt', 'w+')
                    f.write(("EUBED Data for " + self.dosefilename + "\n\n    EUBED for " + r + " = {} " + self.unit + "\n    Mean Dose for " + r + " = {} " + self.unit + "\n    EUBED relative to Mean Dose for " + r + " = {} \n\n").format(EUBED, MEAN, RATIO))
                else:
                    f = open(self.basepath + 'EUBEDData_' + self.dosefilename + '.txt', 'a+')
                    f.write(("    EUBED for " + r + " = {} " + self.unit + "\n    Mean Dose for " + r + " = {} " + self.unit + "\n    EUBED relative to Mean Dose for " + r + " = {} \n\n").format(EUBED, MEAN, RATIO))
            print(("EUBED for " + r + " = {} " + self.patientObject.dcmFileChosen.DoseUnits).format(EUBED))
            print(("Mean Dose for " + r + " = {} " + self.patientObject.dcmFileChosen.DoseUnits).format(MEAN))
            print(("EUBED relative to Mean Dose for " + r + " = {}").format(RATIO))
    
    def EUD(self, ROIList, CreateFile):
        for r in ROIList:
            darr = np.zeros(self.ctObject.quantitiesOfInterest[0].array.shape)
            for i in range(self.ctObject.quantitiesOfInterest[0].array.shape[0]):
                if (i % 30) == 0:
                    prog = i/self.ctObject.quantitiesOfInterest[0].array.shape[0]*100
                    print("Calculating EUD for " + r +  "... (" + str(round(prog,1))+"%)")
                for j in range(self.ctObject.quantitiesOfInterest[0].array.shape[1]):
                    for k in range(self.ctObject.quantitiesOfInterest[0].array.shape[2]):
                        if self.ctObject.structures3D[r][i,j,k] == True:
                            darr[i,j,k] = self.ctObject.quantitiesOfInterest[0].array[i,j,k]
            sum = 0
            if r not in self.TUMORS:
                Alpha_ROI = Alpha_Normal
            elif r in self.TUMORS:
                Alpha_ROI = Alpha_TLiver
            else:
                print('ROI type:' + r + ' not recognized')
                Alpha_ROI = Alpha_Normal
            for i in range(darr.shape[0]):
                for j in range(darr.shape[1]):
                    for k in range(darr.shape[2]):
                        if darr[i,j,k] > 0:
                            sum += math.exp(-Alpha_ROI * darr[i,j,k])
            N = (self.ctObject.structures3D[r]).sum()
            EUD = -1/Alpha_Normal * np.log(sum / N)
            MEAN = np.sum(darr) / N
            RATIO = EUD/MEAN
            if CreateFile == True:
                if ROIList.index(r) == 0:
                    f = open(self.basepath + 'EUDData_' + self.dosefilename + '.txt', 'w+')
                    f.write(("EUD Data for " + self.dosefilename + "\n\n    EUD for " + r + " = {} " + self.unit + "\n    Mean Dose for " + r + " = {} " + self.unit + "\n    EUD relative to Mean Dose for " + r + " = {} \n\n").format(EUD, MEAN, RATIO))
                else:
                    f = open(self.basepath + 'EUDData_' + self.dosefilename + '.txt', 'a+')
                    f.write(("    EUD for " + r + " = {} " + self.unit + "\n    Mean Dose for " + r + " = {} " + self.unit + "\n    EUD relative to Mean Dose for " + r + " = {} \n\n").format(EUD, MEAN, RATIO))
            print(("EUD for " + r + " = {} " + self.patientObject.dcmFileChosen.DoseUnits).format(EUD))
            print(("Mean Dose for " + r + " = {} " + self.patientObject.dcmFileChosen.DoseUnits).format(MEAN))
            print(("EUD relative to Mean Dose for " + r + " = {}").format(RATIO))

    def ConvertDoseUnits(self, seriesdescription = None):
        if self.unit == "Gy/GBq" and str(self.patientObject.dcmFileChosen.DoseUnits) == "Gy/mCi" :
            self.ctObject.quantitiesOfInterest[0].array = 27.027 * self.ctObject.quantitiesOfInterest[0].array
            unitabbr = "GyGBq"
        elif self.unit == "Gy/mCi" and str(self.patientObject.dcmFileChosen.DoseUnits) == "Gy/GBq":
            self.ctObject.quantitiesOfInterest[0].array = self.ctObject.quantitiesOfInterest[0].array / 27.027
            unitabbr = "GymCi"
        else:
            return
        if seriesdescription == None:
            seriesdescription = self.dosefilename + "_" + unitabbr
        name = self.dosefilename + '_' + unitabbr + '.dcm'
        self.ctObject.WriteRTDose(self.ctObject.quantitiesOfInterest[0].array, self.basepath + name, self.unit, seriesdescription)
