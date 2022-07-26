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
from DICOM_RT import DicomPatient as dcmpat
from MIRDCalculation_BED.BioeffectModeling.ROI_Values import *
    
class BioeffectCalculator(dcmpat.PatientCT):
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
        print("ROI's identified:", list(self.ctObject.structures3D.keys()))
        self.BEDimg3D = np.zeros(self.ctObject.img3D.shape)
        

    def BEDCalculator(self):
        if self.unit == "Gy/GBq" and str(self.patientObject.dcmFileChosen.DoseUnits) == "Gy/mCi" :
            self.ctObject.quantitiesOfInterest[0].array = 27.027 * self.ctObject.quantitiesOfInterest[0].array
        elif self.unit == "Gy/mCi" and str(self.patientObject.dcmFileChosen.DoseUnits) == "Gy/GBq":
            self.ctObject.quantitiesOfInterest[0].array = self.ctObject.quantitiesOfInterest[0].array / 27.027
        for i in range(self.ctObject.quantitiesOfInterest[0].array.shape[0]):
            if (i % 20) == 0:
                prog = i/self.ctObject.quantitiesOfInterest[0].array.shape[0]*100
                print("Calculating BED... (" + str(round(prog,1))+"%)")
            else:
                pass
            for j in range(self.ctObject.quantitiesOfInterest[0].array.shape[1]): 
                for k in range(self.ctObject.quantitiesOfInterest[0].array.shape[2]):
                    if self.ctObject.structures3D['Liver'][i,j,k] == True :
                        try:
                            if self.ctObject.structures3D['Right Tumor(s)'][i,j,k] == True or self.ctObject.structures3D['Left Tumor(s)'][i,j,k] == True or self.ctObject.structures3D['All Tumors (Left Lobe)'][i,j,k] == True or self.ctObject.structures3D['All Tumors (Right Lobe)'][i,j,k] == True:
                                Trep = Trep_Tumor
                                AlphaBeta = AlphaBeta_TLiver
                            else:
                                Trep = Trep_Normal
                                AlphaBeta = AlphaBeta_NLiver
                        except:
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

    def EUBEDCalculator(self, ROIList, CreateFile):
        for r in ROIList:
            darr = np.zeros(self.ctObject.quantitiesOfInterest[0].array.shape) 
            for i in range(self.ctObject.quantitiesOfInterest[0].array.shape[0]):
                if (i % 30) == 0:
                    prog = i/self.ctObject.quantitiesOfInterest[0].array.shape[0]*100 
                    print("Calculating EUBED for " + r +  "... (" + str(round(prog,1))+"%)")
                for j in range(self.ctObject.quantitiesOfInterest[0].array.shape[1]): 
                    for k in range(self.ctObject.quantitiesOfInterest[0].array.shape[2]):
                        if self.ctObject.structures3D[r][i,j,k] == True:
                            darr[i,j,k] = self.ctObject.quantitiesOfInterest[0].array[i,j,k]   
            sum = 0
            if r == 'Liver' or r == 'Lung_L' or r == 'Lung_R':
                Alpha_ROI = Alpha_Normal
            elif r == 'Right Tumor(s)' or r == 'Left Tumor(s)' or r == 'All Tumors (Right Lobe)' or r == 'All Tumors (Left Lobe)':
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
            try:
                unit = str(self.patientObject.dcmFileChosen.DoseUnits)
            except:
                unit = 'arb. units'
            if CreateFile == True:
                if ROIList.index(r) == 0:
                    f = open(self.basepath + 'EUBEDData_' + self.dosefilename + '.txt', 'w+')
                    f.write(("EUBED Data for " + self.dosefilename + "\n\n    EUBED for " + r + " = {} " + unit + "\n    Mean Dose for " + r + " = {} " + unit + "\n    EUBED relative to Mean Dose for " + r + " = {} \n\n").format(EUBED, MEAN, RATIO))
                else:
                    f = open(self.basepath + 'EUBEDData_' + self.dosefilename + '.txt', 'a+')
                    f.write(("    EUBED for " + r + " = {} " + unit + "\n    Mean Dose for " + r + " = {} " + unit + "\n    EUBED relative to Mean Dose for " + r + " = {} \n\n").format(EUBED, MEAN, RATIO))
            print(("EUBED for " + r + " = {} " + self.patientObject.dcmFileChosen.DoseUnits).format(EUBED))
            print(("Mean Dose for " + r + " = {} " + self.patientObject.dcmFileChosen.DoseUnits).format(MEAN))
            print(("EUBED relative to Mean Dose for " + r + " = {}").format(RATIO))

class DVH:
    def __init__(self, basepath, dosefile):
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
        print("ROI's identified:", list(self.ctObject.structures3D.keys()))
        self.curves = []
        
    def DVHCalculator(self, ROIList, bins = 1000):
        for r in ROIList:
            darr = np.zeros(self.ctObject.quantitiesOfInterest[0].array.shape) 
            hist = [0] * bins
            maxdose = round(self.ctObject.quantitiesOfInterest[0].array.max())
            mindose = round(self.ctObject.quantitiesOfInterest[0].array.min())
            for i in range(self.ctObject.quantitiesOfInterest[0].array.shape[0]):
                if (i % 30) == 0:
                    prog = i/self.ctObject.quantitiesOfInterest[0].array.shape[0]*100 
                    print("Getting DVH for " + r +  "... (" + str(round(prog,1))+"%)")
                for j in range(self.ctObject.quantitiesOfInterest[0].array.shape[1]): 
                    for k in range(self.ctObject.quantitiesOfInterest[0].array.shape[2]):
                        if self.ctObject.structures3D[r][i,j,k] == True:
                            darr[i,j,k] = self.ctObject.quantitiesOfInterest[0].array[i,j,k]   
            for h in range(bins):
                hist[h] = (darr > (h * (maxdose/bins))).sum() / (darr > 0).sum() * 100
            self.curves.append(hist)
            print(r + ' DVH Calculated.')
    
    def PlotDVHCurves(self, ROIList, bins = 1000):
        if self.ctObject.quantitiesOfInterest[0].array.max() < 0.1:
            print('ERROR: Dose grid values too low to create DVH.')
        try:
            unit = str(self.patientObject.dcmFileChosen.DoseUnits)
        except:
            unit = 'arb. units'
        maxdose = round(self.ctObject.quantitiesOfInterest[0].array.max())
        xpts = np.linspace(0, maxdose, bins)
        for r in ROIList:
            plt.plot(xpts, self.curves[ROIList.index(r)])
        plt.xlabel('Dose [' + unit + ']')
        plt.ylabel("Volume [%]")
        plt.legend(ROIList, loc="upper right")
        plt.savefig(self.basepath + 'DVH_' + self.dosefilename + '.png')
        plt.savefig(self.basepath + 'DVH_' + self.dosefilename + '.jpg')
        plt.show()

