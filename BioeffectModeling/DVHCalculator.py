'''
Plot DVH Curves for selected ROIs

Input args -> DICOM Basepath, RTDOSE filename, List of ROIs to plot, # of bins for DVH

Output -> Matplotlib Plot of DVH curves (saved in png and jpg), DVH points (saved in csv)

'''

import math
import os
import pydicom
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from MIRDCalculation_BED.DICOM_RT import DicomPatient as dcmpat
from MIRDCalculation_BED.BioeffectModeling.ROI_Values import *


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
