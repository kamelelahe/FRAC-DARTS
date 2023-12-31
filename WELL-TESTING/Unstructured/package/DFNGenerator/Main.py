import matplotlib.pyplot as plt
from .Orthogonal import GenerateOrthogonalDFN
from .Conjugate import GenerateConjugateDFN
from .ConjugateSpace import GenerateConjugateDFNWithSpace
from .OrthogonalSpace import GenerateOrthogonalDFNWithSpace
import os

DFN_name = '03Conjugatetest'
# Define output directory
outputDir = 'DFNs/' + str(DFN_name)
# Create directories if they don't exist
os.makedirs(outputDir, exist_ok=True)

type='conjucate'


# Variables for the Generation of the Fracture Networks
Ix=[ 0.0625, 0.125 ]
Iy=[ 0.0625, 0.125 ]
LminX=[6,12]
LminY=[6,12]
proximity_threshold=[0.5]
if type=='orthogonal':
    DFN = GenerateOrthogonalDFNWithSpace(Ix, Iy, LminX, LminY, proximity_threshold=2, outputDir=outputDir, savePic=True)
    #DFN=GenerateOrthogonalDFN(Ix,Iy,LminX,LminY,outputDir, savePic=True)

else:
    DFN= GenerateConjugateDFNWithSpace(Ix, Iy, LminX, LminY,proximity_threshold=proximity_threshold, outputDir=outputDir, savePic=True)
    #DFN = GenerateConjugateDFN(Ix, Iy, LminX, LminY, outputDir, savePic=True)
