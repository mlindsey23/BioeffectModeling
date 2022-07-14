from DICOM_RT.DicomPatient import *

def GetBEDinDICOM(activityObject, ctObject, dicomDirectory, doseFile, structFile):
        if ctObject.img3D.shape[2] < activityObject.img3D.shape[2]:
            activityObject.BEDimg3D = np.zeros((activityObject.img3D.shape[0],activityObject.img3D.shape[1],ctObject.img3D.shape[2]))
        else :
            activityObject.BEDimg3D = np.zeros(activityObject.img3D.shape)
    # LOAD OBJECTS
        dosePath = dicomDirectory + '/' + doseFile
        structPath = dicomDirectory + '/RTSTRUCT/' + structFile
        activityObject.LoadRTDose(dosePath)
        ctObject.LoadStructures(structPath)
    # CALCULATE BED
        for i in range(activityObject.img3D.shape[0]):
            for j in range(activityObject.img3D.shape[1]): 
                if ctObject.img3D.shape[2] < activityObject.img3D.shape[2]:
                    for k in range(ctObject.img3D.shape[2]):
                        if ctObject.structures3D['Liver'][i,j,k] == True :
                            Trep = Trep_Normal
                            AlphaBeta = AlphaBeta_NLiver
                        elif ctObject.structures3D['Lung_L'][i,j,k] == True or ctObject.structures3D['Lung_R'][i,j,k] == True :
                            Trep = Trep_Normal
                            AlphaBeta = AlphaBeta_NLung
                        else :
                            Trep = Trep_Normal
                            AlphaBeta = AlphaBeta_Standard      
                        activityObject.BEDimg3D[i,j,k] = activityObject.img3D[i,j,k] * (1 + (( activityObject.img3D[i,j,k] * Trep) / (AlphaBeta * (Trep + RadionuclideHalfLife))))
                else:
                    for k in range(activityObject.img3D.shape[2]):
                        if ctObject.structures3D['Liver'][i,j,k] == True :
                            Trep = Trep_Normal
                            AlphaBeta = AlphaBeta_NLiver
                        elif ctObject.structures3D['Lung_L'][i,j,k] == True or ctObject.structures3D['Lung_R'][i,j,k] == True :
                            Trep = Trep_Normal
                            AlphaBeta = AlphaBeta_NLung
                        else :
                            Trep = Trep_Normal
                            AlphaBeta = AlphaBeta_Standard
                        activityObject.BEDimg3D[i,j,k] = activityObject.img3D[i,j,k] * (1 + (( activityObject.img3D[i,j,k] * Trep) / (AlphaBeta * (Trep + RadionuclideHalfLife))))
    # PIXEL SCALING
        # Scales max dose = 1
        # Should explore alternative scaling methods
        activityObject.BEDimg3D = activityObject.BEDimg3D / np.max(activityObject.BEDimg3D)
    # RETURN AS DICOM
        ctObject.WriteRTDose(activityObject.BEDimg3D, name = 'BEDCalculation_' + doseFile)
        print('BED Calculated.')
