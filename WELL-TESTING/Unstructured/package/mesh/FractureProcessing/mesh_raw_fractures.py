"""
MIT License
Copyright (c) 2022 Stephan de Hoop 		S.dehoop-1@tudelft.nl
                   Denis Voskov 		D.V.Voskov@tudelft.nl
                   Delft University of Technology, the Netherlands
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from .graph_code import  Graph, create_geo_file
from .calc_intersections_segm_parallel import calc_intersections_segm_parallel
from multiprocessing import Process, freeze_support
import numpy as np
import os
import sys


def mesh_raw_fractures(frac_data_raw, char_len, output_dir='', filename_base='output',
                       height_res=1, apertures_raw=None, box_data=None, margin=25,mesh_raw=False, decimals=5, tolerance_zero=1e-10,
                       tolerance_intersect=1e-10, calc_intersections_before=False, num_partition_x=1,
                       num_partition_y=1, partition_fractures_in_segms=True, matrix_perm=1,
                       char_len_mult=1, char_len_boundary=None, main_algo_iters=2):
    """
    Main fracture preprocessing code, most arguments are optional, but can be tweaked. Please see:
    doi.org/10.1002/essoar.10507519.1 for more explanation on the theory behind the code and some results.
    :param frac_data_raw: input fracture network data in the form [[x1, y1, x2, y2], ..., [...]] [m]
    :param char_len: minimum allowable distance between any two vertices in the domain (characteristic length) [m]
    :param output_dir: directory of the output (cleaned) fracture network and potential .geo and .msh results
    :param filename_base: base name used for all the output files
    :param height_res: height of the resulting 1-layer 3D reservoir [m]
    :param apertures_raw: list of apertures for each fracture segment [m]
    :param box_data: coordinates of the bounding box (in order, from bottom left --> bottom right --> top right --> top left) [m]
    :param margin: margin around the fracture network, used in case no bounding box is given
    :param mesh_raw: boolean, if True then will call GMSH to mesh the raw fracture network in a triangular mesh
    :param decimals: number of decimals used to round-off input data
    :param tolerance_zero: anything below this threshold is considered absolute zero!
    :param calc_intersections_before: boolean, if True calculates intersections between fractures segments before cleaning
    :param num_partition_x: number of partitions used in fracture intersection calculation in x-direction (parallel computing required)
    :param num_partition_y: number of partitions used in fracture intersection calculation in y-direction (parallel computing required)
    :param partition_fractures_in_segms: boolean, if True partitions the fracture network into smaller subsegments of length char_len
    :param matrix_perm: matrix permeability [mD]
    :param char_len_mult: multiplier for mesh characteristic length
    :param char_len_boundary: characteristic mesh length on the boundary (before multiplier)
    :return:
    """
    if apertures_raw is None:
        apertures_raw = np.ones((frac_data_raw.shape[0],)) * 1e-4

    if char_len_boundary is None:
        char_len_boundary = char_len

    if isinstance(frac_data_raw, str):
        frac_data_raw = np.genfromtxt(frac_data_raw)
    assert frac_data_raw.shape[1] == 4, "Data in wrong format, need N rows x 4 columns"

    print('--------------------------------------')
    print('START preprocessing fracture network')
    tot_partitions = num_partition_x * num_partition_y

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    frac_data_raw = np.round(frac_data_raw * 10 ** decimals) * 10 ** (-decimals)
    print('Remove segments of zero length and duplicate segments')
    frac_data_raw, apertures_raw = extract_unique_segms(frac_data_raw, apertures_raw)
    len_raw_sys = np.sqrt((frac_data_raw[:, 0] - frac_data_raw[:, 2]) ** 2 +
                          (frac_data_raw[:, 1] - frac_data_raw[:, 3]) ** 2)

    print('Number of fracture segments: {:}'.format(frac_data_raw.shape[0]))
    print('Min fracture segment length: {:}'.format(np.min(len_raw_sys)))
    print('Max fracture segment length: {:}'.format(np.max(len_raw_sys)))
    print('Mean fracture segment length: {:}'.format(np.mean(len_raw_sys)))
    print('Cleaning length(s): {:}\n'.format(char_len))
    # --------------------------------------------------------------------------

    act_frac_sys = frac_data_raw
    apertures = apertures_raw
    if calc_intersections_before:
        print('START calculating initial intersections raw input fracture network')
        print('\tNOTE: unoptimized!, can take long for very large networks')
        # First find all intersections:
        system_out_par, frac_order_vec_par, partition_lines = \
            calc_intersections_segm_parallel(frac_data_raw, apertures_raw, tolerance_intersect,
                                             tolerance_zero, num_partition_x, num_partition_y)

        # Stack output from all domain together:
        act_frac_sys = system_out_par[0]
        apertures = frac_order_vec_par[0]
        for ii in range(1, tot_partitions):
            act_frac_sys = np.vstack((act_frac_sys, system_out_par[ii]))
            apertures = np.hstack((apertures, frac_order_vec_par[ii]))

        print('DONE calculating initial intersections raw input fracture network')
        if act_frac_sys.shape != frac_data_raw.shape:
            num_intersections = int((act_frac_sys.shape[0] - frac_data_raw.shape[0]) / 2)
            print('\tFound {:} intersections in raw input fracture network\n'.format(num_intersections))
        else:
            print('\tNo intersections found in raw input fracture network\n')

    print('Remove duplicated segments\n')
    act_frac_sys, apertures = extract_unique_segms(act_frac_sys, apertures)
    act_frac_sys_cln = act_frac_sys
    apertures_cln = apertures
    act_frac_sys_raw = np.array(act_frac_sys, copy=True)
    apertures_raw = np.array(apertures_cln, copy=True)

    if np.unique(apertures_cln).size == 1:
        order_cleaning_segms = np.argsort(-np.sqrt((act_frac_sys[:, 0] - act_frac_sys[:, 2])**2 + (act_frac_sys[:, 1] - act_frac_sys[:, 3])**2))
    else:
        order_cleaning_segms = np.argsort(-apertures_cln)

    if partition_fractures_in_segms:
        act_frac_sys_cln, order_cleaning_segms, apertures_cln = segment_fractures(act_frac_sys=act_frac_sys_cln,
                                                                                  char_len=char_len,
                                                                                  order_discr_segms=order_cleaning_segms,
                                                                                  decimals=decimals,
                                                                                  aperture=apertures_cln)

    # --------------------------------------------------------------------------
    print('START constructing graph')
    my_graph = Graph(matrix_perm=matrix_perm)
    my_graph.add_multiple_edges(act_frac_sys_cln)
    my_graph.apertures[np.where(my_graph.active_edges[:my_graph.get_num_edges()] == True)[0]] = apertures_cln
    print('DONE constructing graph\n')

    if box_data is None:
        x_min = np.min(np.min(frac_data_raw[:, [0, 2]])) - margin
        y_min = np.min(np.min(frac_data_raw[:, [1, 3]])) - margin
        x_max = np.max(np.max(frac_data_raw[:, [0, 2]])) + margin
        y_max = np.max(np.max(frac_data_raw[:, [1, 3]])) + margin
        box_data = np.array([[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max]])

    filename_geo_raw = os.path.join(output_dir, filename_base + '_raw_lc_' + str(char_len) + '.geo')
    filename_out_raw = os.path.join(output_dir, filename_base + '_raw_lc_' + str(char_len) + '.msh')
    # --------------------------------------------------------------------------
    print('START writing raw fracture system to file')
    filename_raw = os.path.join(output_dir, filename_base + '_raw_lc_' + str(char_len) + '_fracsys.txt')
    f = open(filename_raw, "w+")
    for frac in act_frac_sys_raw:
        f.write('{:9.5f} {:9.5f} {:9.5f} {:9.5f}\n'.format(frac[0], frac[1], frac[2], frac[3]))
    f.close()

    filename_aper_raw = os.path.join(output_dir, filename_base + '_raw_lc_' + str(char_len) + '_aperture.txt')
    f = open(filename_aper_raw, "w+")
    for aper in apertures_raw:
        f.write('{:16.15f} \n'.format(aper))
    f.close()
    print('DONE writing raw fracture system to file\n')

    print('START creating geo-file for raw network (input for gmsh)')
    create_geo_file(act_frac_sys=act_frac_sys_raw, filename=filename_geo_raw, decimals=decimals,
                    height_res=height_res, char_len=char_len * char_len_mult, box_data=box_data,
                    char_len_boundary=char_len_boundary * char_len_mult)
    print('DONE creating geo-file for raw network (input for gmsh)\n')

    if mesh_raw:
        print('START meshing raw network')
        print(
            '\tNOTE: In gmsh you need to have under Options -> Geometry -> General -> uncheck "Remove duplicate ..." otherwise meshing will crash/take too long')
        os.system("gmsh {:s} -o {:s} -save".format(filename_geo_raw, filename_out_raw))
        print('DONE meshing raw network\n')

    print('Preprocessing succesfully finished')
    print('-----------------------------------')
    return 0


def extract_unique_segms(act_frac_sys, apertures, remove_small_segms=True, tolerance_zero=1e-10):
    """
    Extracts unique fracture segments and apertures
    :param act_frac_sys: input fracture network data in the form [[x1, y1, x2, y2], ..., [...]] [m]
    :param apertures: list of apertures for each fracture segment [m]
    :param remove_small_segms: boolean, if True also removes zero-length fracture segments
    :param tolerance_zero: any value below this will be interpreted as zero!
    :return: unique_segms, unique_apertures
    """
    # Extract unique indices:
    dummy_sys = np.array(act_frac_sys, copy=True)
    indices = np.where(act_frac_sys[:, 0] > act_frac_sys[:, 2])[0]
    dummy_sys[indices, 0:2] = act_frac_sys[indices, 2:]
    dummy_sys[indices, 2:] = act_frac_sys[indices, :2]
    dummy_sys, unique_indices = np.unique(dummy_sys, axis=0, return_index=True)
    dummy_aper = apertures[unique_indices]

    if remove_small_segms:
        len_segm = np.sqrt((dummy_sys[:, 0] - dummy_sys[:, 2])**2 + (dummy_sys[:, 1] - dummy_sys[:, 3])**2)
        dummy_sys = dummy_sys[len_segm > tolerance_zero, :]
        dummy_aper = dummy_aper[len_segm > tolerance_zero]
    return dummy_sys, dummy_aper


def segment_fractures(act_frac_sys, char_len, order_discr_segms, decimals=5, aperture=None):
    """
    Perform partitioning into smaller subsegments (around lc size)
    :param act_frac_sys: input fracture network data in the form [[x1, y1, x2, y2], ..., [...]] [m]
    :param char_len: length of segments after partitioning [m]
    :param order_discr_segms: order in which fracture segments are partitioned and stored in final array
    :param decimals: number of decimals used in rounding off
    :param aperture: array of apertures for each fracture segment
    :return: segmented_fracture_system, new_order_array, new_aperture_array
    """
    length_segms = np.sqrt((act_frac_sys[:, 0] - act_frac_sys[:, 2]) ** 2 +
                           (act_frac_sys[:, 1] - act_frac_sys[:, 3]) ** 2)
    num_new_segms = int(np.sum(np.round(length_segms / char_len)) + np.sum(np.round(length_segms / char_len) == 0))
    act_frac_sys_new = np.zeros((num_new_segms, 4))
    aperture_segm = np.zeros((num_new_segms,))
    order_discr_segms_new = np.zeros((num_new_segms,), dtype=int)
    ith_segm = 0

    for ii in order_discr_segms:
        size_segm = int(max(1, np.round(length_segms[ii] / char_len)))
        id_vec = np.arange(0, size_segm)
        act_frac_sys_new[ith_segm:(ith_segm + size_segm), 0] = act_frac_sys[ii, 0] + id_vec / size_segm * (
                act_frac_sys[ii, 2] - act_frac_sys[ii, 0])
        act_frac_sys_new[ith_segm:(ith_segm + size_segm), 1] = act_frac_sys[ii, 1] + id_vec / size_segm * (
                act_frac_sys[ii, 3] - act_frac_sys[ii, 1])
        act_frac_sys_new[ith_segm:(ith_segm + size_segm), 2] = act_frac_sys[ii, 0] + (id_vec + 1) / size_segm * (
                act_frac_sys[ii, 2] - act_frac_sys[ii, 0])
        act_frac_sys_new[ith_segm:(ith_segm + size_segm), 3] = act_frac_sys[ii, 1] + (id_vec + 1) / size_segm * (
                act_frac_sys[ii, 3] - act_frac_sys[ii, 1])

        if aperture is not None:
            aperture_segm[ith_segm:(ith_segm + size_segm)] = aperture[ii]

        # Update order of fractures:
        order_discr_segms_new[ith_segm:(ith_segm + size_segm)] = np.arange(0, size_segm) + ith_segm

        ith_segm += size_segm

    act_frac_sys_new = np.round(act_frac_sys_new * 10 ** decimals) * 10 ** (-decimals)
    if aperture is not None:
        return act_frac_sys_new, order_discr_segms_new, aperture_segm
    else:
        return act_frac_sys_new, order_discr_segms_new


if __name__ == "__main__":
    freeze_support()
    frac_preprocessing(sys.argv[1], sys.argv[2])