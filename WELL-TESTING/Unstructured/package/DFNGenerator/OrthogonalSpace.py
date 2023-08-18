import numpy as np
from scipy.stats import lognorm
import os
import matplotlib.pyplot as plt
import numba
from numba import jit


@jit(nopython=True)
def check_proximity(frac1, frac2, threshold):
    x1_start, y1_start, x1_end, y1_end = frac1
    x2_start, y2_start, x2_end, y2_end = frac2

    # Calculate the lengths of the fractures
    length_frac1 = ((x1_end - x1_start)**2 + (y1_end - y1_start)**2)**0.5
    length_frac2 = ((x2_end - x2_start)**2 + (y2_end - y2_start)**2)**0.5

    # Determine which fracture is longer
    if length_frac1 > length_frac2:
        long_frac = frac1
        short_frac = frac2
    else:
        long_frac = frac2
        short_frac = frac1

    x1_short=short_frac[0]
    y1_short=short_frac[1]
    x2_short=short_frac[2]
    y2_short=short_frac[3]
    # Segment the longer fracture
    num_segments = int(max(length_frac1, length_frac2) / (threshold * 2))
    if num_segments==0:
        x1_long = long_frac[0]
        y1_long = long_frac[1]
        x2_long= long_frac[2]
        y2_long= long_frac[3]
        condition_1= (x1_long-threshold<=x1_short<=x1_long+threshold) and (y1_long-threshold<=y1_short<=y1_long+threshold)
        condition_2= (x1_long-threshold<=x2_short<=x1_long+threshold) and (y1_long-threshold<=y2_short<=y1_long+threshold)

        condition_3= (x2_long-threshold<=x1_short<=x2_long+threshold) and (y2_long-threshold<=y1_short<=y2_long+threshold)
        condition_4= (x2_long-threshold<=x2_short<=x2_long+threshold) and (y2_long-threshold<=y2_short<=y2_long+threshold)

        # Check if the start or end of the shorter fracture is within the threshold of the current point
        if (condition_1 or condition_2 or condition_3 or condition_4):
            return True
        else:
            return False

    for i in range(num_segments + 1):  # +1 to include the end point
        t = i / num_segments
        x = long_frac[0] + t * (long_frac[2] - long_frac[0])
        y = long_frac[1] + t * (long_frac[3] - long_frac[1])
        condition_1= (x-threshold<=x1_short<=x+threshold) and (y-threshold<=y1_short<=y+threshold)
        condition_2= (x-threshold<=x2_short<=x+threshold) and (y-threshold<=y2_short<=y+threshold)

        # Check if the start or end of the shorter fracture is within the threshold of the current point
        if (condition_1 or condition_2):
            return True

    return False

class GenerateOrthogonalDFNWithSpace:

    def __init__(self, Ix,Iy,LminX,LminY,outputDir,proximity_threshold,savePic=False):
        # -------------------------General Information----------------------------
        cGcell = 2  # matrix cells size
        fGcell = 0.05  # fractures cell size
        cellThickness = 5
        fracThickness = cellThickness

        numCellsX= 250
        numCellsY= 500
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
            maxNumberRealisationsperCase = 1
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
                    horizontalMax=0
                    while p32 <= targetFracIntensity and numIter < maxNumIter and warningCheck < 50000:
                        randInd = np.random.rand()
                        xx = ((1 - randInd) ** (1 / (-2.5 + 1))) * minLength
                        xxx = round(xx / (cGcell + fGcell))
                        while xxx > (numCellsX / 2):
                            randInd = np.random.rand()
                            xx = ((1 - randInd) ** (1 / (-2.5 + 1))) * minLength
                            xxx = round(xx / (cGcell + fGcell))
                        fracLength.append(lengthX[xxx])

                        retry_count = 0
                        max_retries = 100  # Set a limit on the number of retries

                        while retry_count < max_retries:
                            # Randomize y-coordinate
                            posNegGenerator = np.random.randint(0, 2) * 2 - 1
                            yy1 = y1 + round(pd1.rvs()) * posNegGenerator
                            while yy1 > numCellsY / 2 or yy1 < 1:
                                yy1 = y1 + round(pd1.rvs()) * posNegGenerator
                            yCoordFrac1 = np.round(yCoord[yy1 - 1], decimals=4)
                            yCoordFrac2 = yCoordFrac1

                            # Randomize x-coordinate
                            xxGenerator = np.random.randint(0, numCellsX // 2)
                            xCoordFrac1 = np.round(xCoord[xxGenerator] - 10.25, decimals=4)
                            if xCoordFrac1 < minXcoord:
                                xCoordFrac1 = minXcoord
                            xCoordFrac2 = np.round(xCoordFrac1 + fracLength[-1], decimals=4)
                            if xCoordFrac2 > maxXcoord:
                                xCoordFrac2 = maxXcoord
                            if xCoordFrac2 <= xCoordFrac1:
                                xCoordFrac2 = xCoordFrac1 + (fGcell + cGcell)

                            # Check for superimposition and proximity
                            superimpose_condition = any(
                                yCoordFrac1 == existing_frac[1] and (
                                        (xCoordFrac1 >= existing_frac[0] and xCoordFrac1 <= existing_frac[2]) or
                                        (xCoordFrac2 >= existing_frac[0] and xCoordFrac2 <= existing_frac[2])
                                ) for existing_frac in self.horizontal_fractures
                            )
                            proximity_condition = False
                            if not superimpose_condition:
                                proximity_condition = any(
                                    check_proximity((xCoordFrac1, yCoordFrac1, xCoordFrac2, yCoordFrac2), existing_frac,
                                                    proximity_threshold)
                                    for existing_frac in self.horizontal_fractures
                                )

                            if not (superimpose_condition or proximity_condition):
                                break  # Exit the retry loop if conditions are met

                            retry_count += 1

                        if retry_count < max_retries:
                            print('added, Ix= ', str(p32), 'tries= ', str(retry_count))
                            # Corrected fracture length due to violating X coordinates
                            fracLengthCorr = xCoordFrac2 - xCoordFrac1

                            cumArea += fracLengthCorr * fracThickness
                            p32 = cumArea / volRock

                            fracAperture = 1e-4 * (fracLengthCorr ** 0.5)

                            fileID.write(f"{xCoordFrac1} {yCoordFrac1} {xCoordFrac2} {yCoordFrac2}\n")
                            self.horizontal_fractures.append((xCoordFrac1, yCoordFrac1, xCoordFrac2, yCoordFrac2))

                            j += 1
                        else:
                            horizontalMax+=1
                            print( "Max retries reached for horizontal fractures. Ix= ", str(p32),' reached ',str(horizontalMax),' times')


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
                    verticalMax=0
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

                        retry_count = 0
                        max_retries = 100  # Set a limit on the number of retries

                        while retry_count < max_retries:
                            # Randomize x-coordinate
                            posNegGenerator = np.random.randint(0, 2) * 2 - 1
                            xx1 = x1 + round(pd1.rvs()) * posNegGenerator
                            while xx1 > numCellsX / 2 or xx1 < 1:
                                xx1 = x1 + round(pd1.rvs()) * posNegGenerator
                            xCoordFrac1 = np.round(xCoord[xx1 - 1], decimals=4)
                            xCoordFrac2 = xCoordFrac1

                            # Randomize y-coordinate
                            yyGenerator = np.random.randint(0, numCellsY // 2)
                            yCoordFrac1 = np.round(yCoord[yyGenerator] - 10.25, decimals=4)
                            if yCoordFrac1 < minYcoord:
                                yCoordFrac1 = minYcoord
                            yCoordFrac2 = np.round(yCoordFrac1 + fracLength[-1], decimals=4)
                            if yCoordFrac2 > maxYcoord:
                                yCoordFrac2 = maxYcoord
                            if yCoordFrac2 <= yCoordFrac1:
                                yCoordFrac2 = yCoordFrac1 + (fGcell + cGcell)

                            # Check for superimposition and proximity
                            superimpose_condition = any(
                                xCoordFrac1 == existing_frac[0] and (
                                        (yCoordFrac1 >= existing_frac[1] and yCoordFrac1 <= existing_frac[3]) or
                                        (yCoordFrac2 >= existing_frac[1] and yCoordFrac2 <= existing_frac[3])
                                ) for existing_frac in self.vertical_fractures
                            )
                            proximity_condition=False
                            if not superimpose_condition:
                                proximity_condition = any(
                                    check_proximity((xCoordFrac1, yCoordFrac1, xCoordFrac2, yCoordFrac2), existing_frac,
                                                    proximity_threshold)
                                    for existing_frac in self.vertical_fractures
                                )

                            if not (superimpose_condition or proximity_condition):
                                break  # Exit the retry loop if conditions are met

                            retry_count += 1

                        if retry_count < max_retries:
                            print('added, Iy= ', str(p32), 'tries= ', str(retry_count))
                            # Corrected fracture length due to violating Y coordinates
                            fracLengthCorr = yCoordFrac2 - yCoordFrac1

                            cumArea += fracLengthCorr * fracThickness
                            p32 = cumArea / volRock

                            fracAperture = 1e-4 * (fracLengthCorr ** 0.5)

                            fileID.write(f"{xCoordFrac1} {yCoordFrac1} {xCoordFrac2} {yCoordFrac2}\n")
                            self.vertical_fractures.append((xCoordFrac1, yCoordFrac1, xCoordFrac2, yCoordFrac2))
                            j += 1
                        else:
                            verticalMax+=1
                            print( "Max retries reached for vertical fractures. Iy= ", str(p32),' reached ',str(verticalMax),' times')




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
