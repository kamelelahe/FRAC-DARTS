import numpy as np


# Table0-based interpolation procedure
class TableInterpolation:
    def LinearInterP(self, table, x_val, x_index, y_index):
        num = len(table)
        y_val = 0
        for i in range(num):
            if (x_val - table[i][x_index]) <= 0:
                y_val = table[i][y_index] + (x_val - table[i][x_index]) / (table[i][x_index] - table[i - 1][x_index]) \
                        * (table[i][y_index] - table[i - 1][y_index])
                break
        return y_val

    def LinearExtraP(self, table, x_val, x_index, y_index):
        num = len(table)
        y_val = 0
        if (x_val < table[0][x_index]):
            y_val = table[0][y_index] + (x_val - table[0][x_index]) / (table[1][x_index] - table[0][x_index]) * \
                    (table[1][y_index] - table[0][y_index])
        elif (x_val > table[num-1][x_index]):
            y_val = table[num-1][y_index] + (x_val - table[num-1][x_index]) / \
                    (table[num-1][x_index] - table[num-2][x_index]) * (table[num-1][y_index] - table[num-2][y_index])

        return y_val

    def SCALExtraP(self, table, x_val, x_index, y_index):
        num = len(table)
        y_val = 0
        if (x_val < table[0][x_index]):
            y_val = table[0][y_index]
        elif (x_val > table[num - 1][x_index]):
            y_val = table[num - 1][y_index]

        return y_val

        # extrapolate properties if Rs>Rs_max_table
    def SatExtrapolation(self, table, x_val, x_index, y_index, lnum):
        y_val = table[lnum][y_index] + (x_val - table[lnum][x_index]) / \
                (table[lnum][x_index] - table[lnum - 1][x_index]) * (table[lnum][y_index] - table[lnum - 1][y_index])

        return y_val
