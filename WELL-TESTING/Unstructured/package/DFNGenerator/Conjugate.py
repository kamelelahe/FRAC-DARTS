import numpy as np
from scipy.stats import lognorm
import os
import matplotlib.pyplot as plt
import math
class GenerateConjugateDFN:
    def __init__(self, Ix,Iy,LminX,LminY,outputDir, savePic=False):
        # -------------------------General Information----------------------------
        cGcell = 2  # matrix cells size
        fGcell = 0.05  # fractures cell size
        cellThickness = 5
        fracThickness = cellThickness
        """
        numCellsX= 500
        numCellsY= 1000
        minXcoord= 0
        minYcoord= -1024.745
        """
        numCellsX= 500
        numCellsY= 1000
        minXcoord= 0
        minYcoord= 0

        maxXcoord = minXcoord + ((numCellsX / 2) * cGcell) + ((numCellsX / 2) * fGcell)
        maxYcoord = minYcoord + ((numCellsY / 2) * cGcell) + ((numCellsY / 2) * fGcell)

        volRock = (maxXcoord - minXcoord) * (maxYcoord - minYcoord) * cellThickness
        lengthTol = 0.50
        maxNumIter = 10000

        yCoord = [minYcoord + (fGcell / 2) + (i * (cGcell + fGcell)) for i in range(numCellsY // 2)]
        xCoord = [minXcoord + (cGcell + fGcell) * i for i in range(numCellsX // 2)]
        lengthX = [(cGcell + fGcell) * i for i in range(1, numCellsX // 2 + 1)]
        lengthY = [(cGcell + fGcell) * i for i in range(1, numCellsY// 2 + 1)]

        # DoE and related variables
        targetFracIntenDoE = Ix
        targetFracIntenDoEJ = Iy
        LengthMinDoE = LminX
        LengthMinDoEJ = LminY

        # Full Factorial Experimental Design to Sample Combinations of the Fracture Network Variables
        DoE = np.array(np.meshgrid(targetFracIntenDoE, targetFracIntenDoEJ, LengthMinDoE, LengthMinDoEJ)).T.reshape(-1, 4)
        sizeDoE, columnsDoE = DoE.shape
        number = 0

        for DoELoop in range(sizeDoE):

            targetFracIntensity = DoE[DoELoop][0]
            targetFracIntensityJ = DoE[DoELoop][1]
            minLength = DoE[DoELoop][2]
            minLengthJ = DoE[DoELoop][3]
            maxNumberRealisationsperCase = 2
            self.horizontal_fractures = []  # Each element is a tuple: (x1, y1, x2, y2)
            self.vertical_fractures = []  # Each element is a tuple: (x1, y1, x2, y2)

            for multRealisations in range(maxNumberRealisationsperCase):

                PropertyDir= outputDir + '/properties'
                os.makedirs(PropertyDir, exist_ok=True)
                propertiesFile=f"{PropertyDir}\\DFNProperties{number + 1:03}.txt"
                with open(propertiesFile, 'w') as fileID:
                    fileID.write(f"Ix={targetFracIntensity}\n"
                                 f"Iy= {targetFracIntensityJ}\n"
                                 f"minLengthX {minLength}\n"
                                 f"minLengthY {minLengthJ}\n"
                                 f"Realizarion {multRealisations}\n")

                TextDir= outputDir + '/Text'
                os.makedirs(TextDir, exist_ok=True)
                filename = f"{TextDir}\\DFN{number + 1:03}.txt"
                with open(filename, 'w') as fileID:

                    y1 = np.random.randint(20, numCellsY // 2 - 20)

                    mean1 = 125
                    SD1 = mean1 * 1.6
                    meanLog1 = np.log(mean1 ** 2 / np.sqrt(SD1 ** 2 + mean1 ** 2))
                    SDLog1 = np.sqrt(np.log((SD1 ** 2 / mean1 ** 2) + 1))
                    pd1 = lognorm(s=SDLog1, scale=np.exp(meanLog1))

                    cumArea = 0
                    p32 = 0
                    j = 0  # Python uses 0-based indexing
                    m = len(yCoord)
                    numberHits = np.zeros(m)
                    counter1 = 0
                    numIter = 0
                    warningCheck = 0

                    fracLength = []

                    while p32 <= targetFracIntensity and numIter < maxNumIter and warningCheck < 50000:
                        randInd = np.random.rand()
                        xx = ((1 - randInd) ** (1 / (-2.5 + 1))) * minLength
                        xxx = round(xx / (cGcell + fGcell))
                        while xxx > (numCellsX / 2):
                            randInd = np.random.rand()
                            xx = ((1 - randInd) ** (1 / (-2.5 + 1))) * minLength
                            xxx = round(xx / (cGcell + fGcell))
                        fracLength.append(lengthX[xxx])

                        posNegGenerator = np.random.randint(0, 2) * 2 - 1
                        yy1 = y1 + round(pd1.rvs()) * posNegGenerator
                        warningCheck = 0
                        while yy1 > numCellsY / 2 or yy1 < 1:
                            yy1 = y1 + round(pd1.rvs()) * posNegGenerator
                            warningCheck += 1
                        if warningCheck > 49999:
                            print('problems with sampling in X')
                        yCoordFrac1 = np.round(yCoord[yy1 - 1], decimals=4)

                        xxGenerator = np.random.randint(0, numCellsX // 2)
                        xCoordFrac1 = np.round(xCoord[xxGenerator] - 10.25, decimals=4)

                        # Calculate the end points based on the 30-degree angle
                        xCoordFrac2 = np.round(xCoordFrac1 + fracLength[-1] * np.cos(np.radians(30)), decimals=4)
                        yCoordFrac2 = np.round(yCoordFrac1 + fracLength[-1] * np.sin(np.radians(30)), decimals=4)

                        # Check if the coordinates are within the domain
                        if xCoordFrac2 > maxXcoord or yCoordFrac2 > maxYcoord:
                            continue

                        # Check for superimposition if there's more than one fracture
                        if j > 0:
                            loopCounter = 0
                            while loopCounter < j and numIter < maxNumIter:
                                # Check if fractures superimpose in space
                                superimpose_condition = (
                                        yCoordFrac1 == self.horizontal_fractures[loopCounter][1]
                                        and (
                                                (xCoordFrac1 >= self.horizontal_fractures[loopCounter][
                                                    0] and xCoordFrac1 <= self.horizontal_fractures[loopCounter][2])
                                                or (xCoordFrac2 >= self.horizontal_fractures[loopCounter][
                                            0] and xCoordFrac2 <= self.horizontal_fractures[loopCounter][2])
                                        )
                                )
                                if superimpose_condition:
                                    # Resample the x-coordinates
                                    xxGenerator = np.random.randint(0, numCellsX // 2)
                                    xCoordFrac1 = np.round(xCoord[xxGenerator] - 10.25, decimals=4)
                                    xCoordFrac2 = np.round(xCoordFrac1 + fracLength[-1] * np.cos(np.radians(30)),
                                                           decimals=4)
                                    yCoordFrac2 = np.round(yCoordFrac1 + fracLength[-1] * np.sin(np.radians(30)),
                                                           decimals=4)

                                    # Reset the loopCounter to check again
                                    loopCounter = 0
                                    numIter += 1
                                else:
                                    loopCounter += 1

                        # Corrected fracture length due to violating X coordinates
                        fracLengthCorr = np.sqrt((xCoordFrac2 - xCoordFrac1) ** 2 + (yCoordFrac2 - yCoordFrac1) ** 2)

                        cumArea += fracLengthCorr * fracThickness
                        p32 = cumArea / volRock

                        fracAperture = 1e-4 * (fracLengthCorr ** 0.5)

                        fileID.write(f"{xCoordFrac1} {yCoordFrac1} {xCoordFrac2} {yCoordFrac2}\n")
                        self.horizontal_fractures.append((xCoordFrac1, yCoordFrac1, xCoordFrac2, yCoordFrac2))

                        j += 1
                    x1 = np.random.randint(20, numCellsX // 2 - 20)

                    mean1 = 125
                    SD1 = mean1 * 1.6
                    meanLog1 = np.log(mean1 ** 2 / np.sqrt(SD1 ** 2 + mean1 ** 2))
                    SDLog1 = np.sqrt(np.log((SD1 ** 2 / mean1 ** 2) + 1))
                    pd1 = lognorm(s=SDLog1, scale=np.exp(meanLog1))

                    cumArea = 0
                    p32 = 0
                    j = 0  # Python uses 0-based indexing
                    m = len(xCoord)
                    numberHits = np.zeros(m)
                    counter1 = 0
                    numIter = 0
                    warningCheck = 0

                    fracLength = []

                    while p32 <= targetFracIntensityJ and numIter < maxNumIter and warningCheck < 50000:
                        randInd = np.random.rand()
                        yy = ((1 - randInd) ** (1 / (-2.5 + 1))) * minLengthJ
                        yyy = round(yy / (cGcell + fGcell))
                        while yyy > (numCellsY / 2):
                            randInd = np.random.rand()
                            yy = ((1 - randInd) ** (1 / (-2.5 + 1))) * minLengthJ
                            yyy = round(yy / (cGcell + fGcell))
                        fracLength.append(lengthY[yyy])

                        posNegGenerator = np.random.randint(0, 2) * 2 - 1

                        xx1 = x1 + round(pd1.rvs()) * posNegGenerator
                        warningCheck = 0
                        while xx1 > numCellsX / 2 or xx1 < 1:
                            xx1 = x1 + round(pd1.rvs()) * posNegGenerator
                            warningCheck += 1
                        if warningCheck > 49999:
                            print('problems with sampling in Y')
                        xCoordFrac1 = np.round(xCoord[xx1 - 1], decimals=4)
                        xCoordFrac2 = xCoordFrac1

                        yyGenerator = np.random.randint(0, numCellsY // 2)
                        yCoordFrac1 = np.round(yCoord[yyGenerator] - 10.25, decimals=4)

                        if yCoordFrac1 < minYcoord:
                            yCoordFrac1 = minYcoord

                        yCoordFrac2 = np.round(yCoordFrac1 + fracLength[-1], decimals=4)

                        if yCoordFrac2 > maxYcoord:
                            yCoordFrac2 = maxYcoord

                        if yCoordFrac2 <= yCoordFrac1:
                            yCoordFrac2 = yCoordFrac1 + (fGcell + cGcell)
                        # Check for superimposition if there's more than one fracture
                        if j > 0:
                            loopCounter = 0
                            while loopCounter < j and numIter < maxNumIter:
                                # Check if fractures superimpose in space
                                superimpose_condition = (
                                        xCoordFrac1 == self.vertical_fractures[loopCounter][0]
                                        and (
                                                (yCoordFrac1 >= self.vertical_fractures[loopCounter][
                                                    1] and yCoordFrac1 <= self.vertical_fractures[loopCounter][3])
                                                or (yCoordFrac2 >= self.vertical_fractures[loopCounter][
                                            1] and yCoordFrac2 <= self.vertical_fractures[loopCounter][3])
                                        )
                                )
                                if superimpose_condition:
                                    # Resample the y-coordinates
                                    yyGenerator = np.random.randint(0, numCellsY // 2)
                                    yCoordFrac1 = np.round(yCoord[yyGenerator] - 10.25, decimals=4)
                                    if yCoordFrac1 < minYcoord:
                                        yCoordFrac1 = minYcoord
                                    yCoordFrac2 = np.round(yCoordFrac1 + fracLength[-1], decimals=4)
                                    if yCoordFrac2 > maxYcoord:
                                        yCoordFrac2 = maxYcoord
                                    if yCoordFrac2 <= yCoordFrac1:
                                        yCoordFrac2 = yCoordFrac1 + (fGcell + cGcell)

                                    # Reset the loopCounter to check again
                                    loopCounter = 0
                                    numIter += 1
                                else:
                                    loopCounter += 1

                        # Corrected fracture length due to violating Y coordinates
                        fracLengthCorr = yCoordFrac2 - yCoordFrac1

                        cumArea += fracLengthCorr * fracThickness
                        p32 = cumArea / volRock

                        fracAperture = 1e-4 * (fracLengthCorr ** 0.5)

                        fileID.write(f"{xCoordFrac1} {yCoordFrac1} {xCoordFrac2} {yCoordFrac2}\n")
                        self.vertical_fractures.append((xCoordFrac1, yCoordFrac1, xCoordFrac2, yCoordFrac2))
                        j += 1

            if savePic:
                ImgDir= outputDir + '/Images'
                os.makedirs(ImgDir, exist_ok=True)
                # Plotting
                plt.figure(figsize=(7, 14))

                # Plot horizontal fractures in blue
                for frac in self.horizontal_fractures:
                    plt.plot([frac[0], frac[2]], [frac[1], frac[3]], color='blue', linewidth=0.7)

                # Plot vertical fractures in red
                for frac in self.vertical_fractures:
                    plt.plot([frac[0], frac[2]], [frac[1], frac[3]], color='red', linewidth=0.7)

                plt.xlabel('X Coordinate')
                plt.ylabel('Y Coordinate')
                plt.title('Fracture Visualization')
                plt.grid(True)
                imagename = f"{ImgDir}\\DFN{number + 1:03}.png"
                plt.savefig(imagename, format='png', dpi=300)

            number+=1
