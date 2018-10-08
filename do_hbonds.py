#!/usr/bin/env python

# This script automates hbonds analysis with gromacs. By interpreting .xpm matrix
# provides 2 concise .txt summaries and lots of .csv tables

from __future__ import print_function, division
import sys
import numpy as np
import os
import argparse
import uuid
import subprocess
import re
import shlex


def read_xvg(file_handle):
    """Parses XVG file legends and data"""

    _ignored = {'legend', 'view'}
    _re_series = re.compile('s[0-9]+$')
    _re_xyaxis = re.compile('[xy]axis$')

    metadata = {}
    num_data = []

    metadata['labels'] = {}
    metadata['labels']['series'] = []

    for line in file_handle:
        line = line.strip()
        if line.startswith('@'):
            tokens = shlex.split(line[1:])
            if tokens[0] in _ignored:
                continue
            elif tokens[0] == 'TYPE':
                if tokens[1] != 'xy':
                    raise ValueError('Chart type unsupported: \'{0}\'. Must be \'xy\''.format(tokens[1]))
            elif _re_series.match(tokens[0]):
                metadata['labels']['series'].append(tokens[-1])
            elif _re_xyaxis.match(tokens[0]):
                metadata['labels'][tokens[0]] = tokens[-1]
            elif len(tokens) == 2:
                metadata[tokens[0]] = tokens[1]
            else:
                print('Unsupported entry: {0} - ignoring'.format(tokens[0]), file=sys.stderr)
        elif line[0].isdigit():
                num_data.append(map(float, line.split()))

    num_data = list(zip(*num_data))

    if not metadata['labels']['series']:
        for series in range(len(num_data) - 1):
            metadata['labels']['series'].append('')

    return metadata, num_data


def unquote(s):
    return s[1+s.find('"'):s.rfind('"')]


def uncomment(s):
    return s[3+s.find('/*'):s.rfind('*/')-1]


def make_lut(c, nb):
    color = c.split('/*')
    value = unquote(color[1])
    if value == "None":
        value = 0.0
    if value == "Present":
        value = 1.0
    tmp = unquote(color[0]).split(' c ')
    key = c[1:1+nb]
    color = tmp[1].strip()
    return key, (color, float(value))


def process_meta(metadata):
    uncommented = []
    for i in metadata:
        j = uncomment(i)
        if len(j) == 0:
            continue
        if j[0] == ' ':
            uncommented[len(uncommented)-1] += j
        elif ":" in j:
            uncommented.append(j)
    m_dict = dict()
    for i in uncommented:
        index = i.find(':')
        key = i[:index]
        value = i[index+1:].strip()
        if value[0] == '"':
            value = unquote(value)
        if key in m_dict:
            m_dict[key] += ' '+value
        else:
            m_dict[key] = value
    return m_dict


def read_xpm(fhandle, reverse=False):

        # Read in lines until we fidn the start of the array
        metadata = [fhandle.readline()]
        while not metadata[-1].startswith("static char *gromacs_xpm[]"):
            metadata.append(fhandle.readline())

        # The next line will contain the dimensions of the array
        dim = fhandle.readline()
        # There are four integers surrounded by quotes
        nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]

        # The next dim[2] lines contain the color definitions
        # Each pixel is encoded by dim[3] bytes, and a comment
        # at the end of the line contains the corresponding value
        lut = dict([make_lut(fhandle.readline(), nb) for _ in range(nc)])

        colors = []
        values = []

        for i in fhandle:
            if i.startswith("/*"):
                metadata.append(i)
                continue
            j = unquote(i)
            colors.append([lut[j[k:k + nb]][0] for k in range(0, nx, nb)])
            values.append([lut[j[k:k + nb]][1] for k in range(0, nx, nb)])
        metadata = process_meta(metadata)

        if reverse:
            values.reverse()
            colors.reverse()

        return colors, values, metadata, lut


def read_index(fhandle):
    """Reads an already open index file and returns a list of tuples (groupname, list fo lines)"""
    result = []
    current_label = None
    current_value = []

    for line in fhandle:
        line = line.strip()
        if line.startswith('[') and line.endswith(']'):
            if current_label:
                result.append((current_label, current_value))
                current_value = []
            current_label = line[1:-1]
            continue
        if len(line) == 0:
            continue
        current_value.append(line)
    if current_label:
        result.append((current_label, current_value))
    else:
        raise IOError('File {0} not readble: missing [group]'.format(fhandle.name))
    return result


def read_gro(fhandle):
    """Reads and already open .gro file and returns a tuple of 0 indexed arrays resname resnumber and atomname"""
    gro_lines = fhandle.readlines()
    gro_lines.pop(0)  # remove comment
    n_atoms = int(gro_lines.pop(0))
    gro_lines.pop(-1)  # remove last line containing the box
    # box = [float(i) for i in box_list]
    gro_atoms = dict()
    gro_residues = dict()
    res_adder = 0
    atom_adder = 0
    res_step = 100000
    atom_step = 100000
    for line in gro_lines:
        r_nam = line[5:10].strip()
        rn = line[0:5]
        r_num = str(int(rn) + res_adder)
        a_nam = line[10:15].strip()
        an = line[15:20]
        a_num = str(int(an) + atom_adder)
        if rn == '99999':
            res_adder += res_step
        if an == '99999':
            atom_adder += atom_step
        gro_atoms[a_num] = {'aname': a_nam, 'rname': r_nam, 'rnum': r_num}
        if r_num in gro_residues:
            gro_residues[r_num]['atoms'].append(a_num)
        else:
            gro_residues[r_num] = {'rname': r_nam, 'atoms': [a_num]}
    assert n_atoms == len(gro_atoms)
    return gro_atoms, gro_residues


def backup_if_exists(file_name):
    if os.path.exists(file_name):
        file_name = os.path.abspath(file_name)
        the_dir = os.path.dirname(file_name)
        base = os.path.basename(file_name)
        files = os.listdir(the_dir)
        backups = [i for i in files if i[:len(base)+1] == ('#' + base)]
        if len(backups) == 0:
            number = 1
        else:
            number = max([int(i.split('.')[-1][:-1]) for i in backups])+1
        new_name = os.path.join(the_dir, '#'+base+'.'+str(number)+'#')
        os.rename(file_name, new_name)
        file_name = os.path.join(the_dir, os.path.basename(file_name))
        new_name = os.path.join(the_dir, os.path.basename(new_name))
        sys.stdout.write('\nBack Off! I just backed up {} to {}\n\n'.format(file_name, new_name))


def delete_if_exists(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)


def check_read_ext(choices):
    class Act(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            extension = os.path.splitext(values)[1][1:]
            c = choices
            if type(c) == str:
                c = {c}
            if extension not in c:
                parser.error("file {} doesn't end with one of {}".format(values, choices))
            elif os.path.exists(values):
                if os.path.isdir(values):
                    parser.error("{} is a directory".format(values))
                elif not os.access(values, os.R_OK):
                    parser.error("file {} is not readable".format(values))
                else:
                    setattr(namespace, self.dest, values)
            else:
                parser.error("file {} doesn't exists".format(values))
    return Act


def check_write_ext(choices):
    class Act(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            extension = os.path.splitext(values)[1][1:]
            c = choices
            if type(c) == str:
                c = {c}
            if extension not in c:
                parser.error("file {} doesn't end with one of {}".format(values, choices))
            elif os.path.exists(values):
                if os.path.isdir(values):
                    parser.error("{} is a directory".format(values))
                elif not os.access(values, os.W_OK):
                    parser.error("file {} is not writeable".format(values))
                else:
                    setattr(namespace, self.dest, values)
            else:
                path = os.path.dirname(os.path.abspath(values))
                if not os.access(path, os.W_OK):
                    parser.error("the directory {} is not writeable".format(path))
                else:
                    setattr(namespace, self.dest, values)
    return Act


def check_write_dir():
    class Act(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if os.path.exists(values):
                if os.path.isfile(values):
                    parser.error("{} is a file".format(values))
                elif not os.access(values, os.W_OK):
                    parser.error("directory {} is not writeable".format(values))
                else:
                    setattr(namespace, self.dest, values)
            else:
                path = os.path.dirname(os.path.abspath(values))
                if not os.access(path, os.W_OK):
                    parser.error("the directory {} is not writeable".format(path))
                else:
                    setattr(namespace, self.dest, values)
    return Act


if __name__ == '__main__':

    ######################
    #                    #
    # Parse command line #
    #                    #
    ######################
    p = argparse.ArgumentParser()
    p.add_argument('-f', action=check_read_ext({'xtc', 'trr'}), required=True, dest='trajectory_file',
                   help='trajectory file .trr or .xtc')
    p.add_argument('-s', action=check_read_ext('tpr'), required=True, dest='topology_file',
                   help='topology file .tpr')
    p.add_argument('-o', action=check_write_dir(), required=True, dest='output_dir',
                   help='output directory')
    p.add_argument('-sel', type=int, required=True, dest='selection',
                   help='selection group 1')
    # p.add_argument('-g', action=check_write_ext('ndx'), required=False, dest='group_index',
    #               help='keep the group index file, defaults to hbond_group.ndx', default='hbond_group.ndx')
    p.add_argument('-sel2', type=int, required=False, dest='selection2',
                   help='selection group 2, defaults to selection group 1')
    p.add_argument('-l', nargs='?', action=check_write_ext('txt'), required=False, dest='list_file',
                   help='hbond list file, defaults to hbond_list.txt', default='hbond_list.txt')
    p.add_argument('-n', action=check_read_ext('ndx'), required=False, dest='index_file',
                   help='index file')
    p.add_argument('-dt', type=int, required=False,
                   help=' selection2')
    p.add_argument('-b', nargs='?', type=int, dest='begin', required=False, default=0,
                   help='begin time for existence check (ps), defaults to 0')
    p.add_argument('-e', type=int, dest='end', required=False, help='end time for existence check (ps)')
    p.add_argument('-t', nargs='?', type=int, dest='threshold', required=False,
                   help='occupancy threshold to accept hbond (percent), defaults to 0', default=0)
    p.add_argument('-min-res-dist', nargs='?', dest='min_res_dist', required=False, type=int, default=0,
                   help='consider hbonds with residue number difference greater than, defaults to 0 ')
    args = p.parse_args()

    ###############################################
    #                                             #
    # Run gmx hbond to obtain .ndx and .xpm files #
    #                                             #
    ###############################################
    trajectory_file = os.path.abspath(args.trajectory_file)
    topology_file = os.path.abspath(args.topology_file)
    output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    index_string = ''
    # generate random names for matrix and hbond index files
    # matrix_file = os.path.join(output_dir, str(uuid.uuid4())+'.xpm')
    # hbond_index_file = os.path.join(output_dir, str(uuid.uuid4())+'.ndx')
    matrix_file = os.path.join(output_dir, 'hbmap.xpm')
    hbond_index_file = os.path.join(output_dir, 'hbond.ndx')
    number_file = os.path.join(output_dir, str(uuid.uuid4())+'.xvg')
    sel1 = str(args.selection)
    sel2 = str(args.selection)
    if args.selection2:
        sel2 = str(args.selection2)
    command = ['gmx', 'hbond', '-f', trajectory_file, '-s', topology_file, '-hbn', hbond_index_file,
               '-hbm', matrix_file, '-num', number_file]
    if args.index_file:
        command.append('-n')
        command.append(os.path.abspath(args.index_file))
    if args.dt:
        command.append('-dt')
        command.append(str(args.dt))
    input_string = sel1 + '\n' + sel2 + '\n'
    sys.stdout.write('\nStarting gmx hbond\n\nRunning, please wait...\n\n')
    gmx = subprocess.run(command, input=input_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    delete_if_exists(number_file)
    if gmx.returncode:
        delete_if_exists(matrix_file)
        delete_if_exists(hbond_index_file)
        sys.stderr.write('\nERROR: called Gromacs with line:\n{}\n\n'.format(' '.join(command)))
        sys.stderr.write('INPUT:\n{}\n\n'.format(input_string))
        sys.stderr.write('RETURN CODE: {}\n\n'.format(gmx.returncode))
        sys.stderr.write('STDOUT:\n{}\n\n'.format(gmx.stdout))
        sys.stderr.write('STERR:\n{}\n\n'.format(gmx.stderr))
        sys.exit(1)
    sys.stdout.write('{}\n'.format(gmx.stdout))
    sys.stdout.write('Successfully completed gmx hbond\n\n')

    ######################################
    #                                    #
    # generate the .gro from gmx trjconv #
    #                                    #
    ######################################
    gro_file = os.path.join(output_dir, str(uuid.uuid4())+'.gro')
    command = ['gmx', 'trjconv', '-f', trajectory_file, '-s', topology_file, '-b', '0', '-e', '0', '-o', gro_file]
    input_string = '0\n'
    sys.stdout.write('Starting gmx trjconv\n\n')
    gmx = subprocess.run(command, input=input_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    if gmx.returncode:
        delete_if_exists(matrix_file)
        delete_if_exists(hbond_index_file)
        delete_if_exists(gro_file)
        sys.stderr.write('\nERROR: called Gromacs with line:\n{}\n\n'.format(' '.join(command)))
        sys.stderr.write('INPUT:\n{}\n\n'.format(input_string))
        sys.stderr.write('RETURN CODE: {}\n\n'.format(gmx.returncode))
        sys.stderr.write('STDOUT:\n{}\n\n'.format(gmx.stdout))
        sys.stderr.write('STERR:\n{}\n\n'.format(gmx.stderr))
        sys.exit(1)
    sys.stdout.write('{}\n'.format(gmx.stdout))
    sys.stdout.write('Successfully completed gmx trjconv\n\n')

    ##################################################
    #                                                #
    # read .gro, .ndx and .xpm files and delete them #
    #                                                #
    ##################################################
    with open(gro_file, 'r') as gro_handle:
        atoms, _ = read_gro(gro_handle)
    with open(hbond_index_file) as index_handle:
        indexes = read_index(index_handle)
    with open(matrix_file) as matrix_handle:
        _, matrix, meta, _ = read_xpm(matrix_handle, reverse=True)
    os.remove(gro_file)
    # os.remove(hbond_index_file)
    # os.remove(matrix_file)

    ###########################
    #                         #
    # perform hbond selection #
    #                         #
    ###########################

    # select relevant columns of the xpm matrix
    x_values = meta['x-axis'].split()
    x_values = [int(i) for i in x_values]
    if not args.end:
        end = x_values[-1]
    else:
        end = int(args.end)
    start_time = 0
    end_time = len(x_values) - 1
    for bond_number in range(len(x_values)):
        if x_values[bond_number] >= args.begin:
            start_time = bond_number
            break
    for bond_number in reversed(range(len(x_values))):
        if x_values[bond_number] <= end:
            end_time = bond_number
            break

    # atom positions are the last group of the index
    _, lines = indexes[-1]

    # we want columns along time
    consensus = np.array(matrix).transpose()

    occupancy = np.sum(consensus[start_time:end_time, :], axis=0) / (end_time - start_time)
    threshold = float(args.threshold) / 100.0

    # open list file #
    if len(os.path.dirname(args.list_file)) != 0:
        list_file = os.path.abspath(args.list_file)
    else:
        list_file = os.path.join(output_dir, args.list_file)
    backup_if_exists(list_file)
    list_handle = open(list_file, "w")
    list_handle.write('; #  |          Donor          |  Hydrogen   |         Acceptor        | occup.\n')
    list_handle.write('-------------------------------------------------------------------------------\n')
    # open groups file
    # if len(os.path.dirname(args.group_index)) != 0:
    #     groups_file = os.path.abspath(args.group_index)
    # else:
    #     groups_file = os.path.join(output_dir, args.group_index)
    groups_file = os.path.join(output_dir, 'group_index.ndx')
    backup_if_exists(groups_file)
    groups_handle = open(groups_file, 'w')

    select_string = ''
    columns = []
    labels = []
    bond_index = 0
    bond_number = 0
    for the_line, the_occupancy in zip(lines, occupancy):
        # TODO remove if true
        if True:  # the_occupancy > threshold:
            dnum, hnum, anum = the_line.strip().split()
            dname = atoms[dnum]['aname']
            drname = atoms[dnum]['rname']
            drnum = atoms[dnum]['rnum']
            hname = atoms[hnum]['aname']
            aname = atoms[anum]['aname']
            arname = atoms[anum]['rname']
            arnum = atoms[anum]['rnum']
            if abs(int(drnum) - int(arnum)) > args.min_res_dist and the_occupancy > threshold:
                list_handle.write("{:>4} | ".format(str(bond_number)))
                list_handle.write("{:>4} {:6} {:>4} {:6} | ".format(drname, drnum, dname, dnum))
                list_handle.write("{:>4} {:6} | ".format(hname, hnum))
                list_handle.write("{:>4} {:6} {:>4} {:6} | {:>3}%\n".format(arname, arnum, aname, anum,
                                                                            int(round(the_occupancy * 100.0, 0))))
                columns.append(bond_index)
                labels.append("{}{}-{}{}".format(drname, drnum, arname, arnum))
                # groups_handle.write("[ {}{}_{}{} ]\n".format(drname, drnum, arname, arnum))
                # groups_handle.write("{0:>6} {1:>6}\n".format(dnum, anum))
                # select_string += str(bond_number) + '\n'
                bond_number += 1
            select_string += str(bond_index) + '\n'
            groups_handle.write("[ {}{}_{}{} ]\n".format(drname, drnum, arname, arnum))
            groups_handle.write("{0:>6} {1:>6}\n".format(dnum, anum))
        bond_index += 1
    groups_handle.close()
    list_handle.close()
    sys.stdout.write('File {} written\n\n'.format(groups_file))
    sys.stdout.write('File {} written\n\n'.format(list_file))

    ####################
    #                  #
    # run gmx distance #
    #                  #
    ####################
    # generate random name for distance file
    # distance_file = os.path.join(output_dir, str(uuid.uuid4()) + '.xvg')
    distance_file = os.path.join(output_dir, 'hbond_dist.xvg')
    command = ['gmx', 'distance', '-f', trajectory_file, '-s', topology_file, '-oav', distance_file, '-n',
               groups_file]
    if args.dt:
        command.append('-dt')
        command.append(str(args.dt))
    sys.stdout.write('Starting gmx distance\n\nRunning, please wait...\n\n')
    gmx = subprocess.run(command, input=select_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True)
    if gmx.returncode:
        delete_if_exists(distance_file)
        sys.stderr.write('\nERROR: called Gromacs with line:\n{}\n\n'.format(' '.join(command)))
        sys.stderr.write('RETURN CODE: {}\n\n'.format(gmx.returncode))
        sys.stderr.write('STDOUT:\n{}\n\n'.format(gmx.stdout))
        sys.stderr.write('STERR:\n{}\n\n'.format(gmx.stderr))
        sys.exit(1)
    # sys.stdout.write('{}\n'.format(gmx.stdout))
    sys.stdout.write('Successfully completed gmx distance\n\n')

    ####################################
    #                                  #
    # read distance file and delete it #
    #                                  #
    ####################################
    with open(distance_file, 'r') as distance_handle:
        _, distance_data = read_xvg(distance_handle)
    # os.remove(distance_file)

    # TODO refactor read_xvg to avoid this transpose
    distances = np.array(distance_data).transpose()

    time = distances[:, 0]
    distances = distances[:, 1:]

    # assert(distances.shape == consensus.shape)

    consensus = consensus[:, columns]
    distances = distances[:, columns]

    # assert(distances.shape == consensus.shape)

    ##########################################
    #                                        #
    # collapse hbonds between same residuals #
    #                                        #
    ##########################################
    non_unique_pairs = []  # list of [label, start_column, extent]
    for i, label in enumerate(labels):
        if len(non_unique_pairs) == 0:
            non_unique_pairs.append([label, i, 1])
        elif non_unique_pairs[-1][0] == label:
            non_unique_pairs[-1][2] += 1
        else:
            non_unique_pairs.append([label, i, 1])

    unique_pairs = [non_unique_pairs.pop(0)]  # list of [label, startcolumn1, extent1, startcolumn2, extent2, ...]
    while len(non_unique_pairs) > 0:
        pair1 = non_unique_pairs.pop(0)
        found = False
        for i, pair2 in enumerate(unique_pairs):
            if pair1[0] == pair2[0]:
                unique_pairs[i].append(pair1[1])
                unique_pairs[i].append(pair1[2])
                found = True
                break
        if not found:
            unique_pairs.append(pair1)

    nrow = distances.shape[0]
    ncol = len(unique_pairs)

    distances1 = np.zeros((nrow, ncol), dtype='float32')
    consensus1 = np.zeros((nrow, ncol), dtype='float32')

    for col, pair in enumerate(unique_pairs):
        n_chunks = (len(pair) - 1) // 2
        for row in range(nrow):
            dist = 1.0e20
            cons = 0.0
            for i in range(n_chunks):
                start_column = pair[1 + i * 2]
                extent = pair[2 + i * 2]
                d1 = distances[row, start_column:start_column + extent].min()
                dist = min(dist, d1)
                cons = cons + consensus[row, start_column:start_column + extent].sum()
            distances1[row, col] = dist
            if cons > 0.0:
                consensus1[row, col] = 1.0
            else:
                consensus1[row, col] = 0.0

    #######################
    #                     #
    # write all the files #
    #                     #
    #######################
    for col, pair in enumerate(unique_pairs):
        bond = pair[0]
        f_name = os.path.join(output_dir, 'hbond-' + bond + '.csv')
        backup_if_exists(f_name)
        with open(f_name, 'w') as f_handle:
            f_handle.write('Time,Distance,Consensus\n')
            for row in range(nrow):
                if consensus1[row, col]:
                    cons = 'Present'
                else:
                    cons = 'None'
                f_handle.write('{0:},{1:.3f},{2:}\n'.format(time[row], distances1[row, col], cons))
        sys.stdout.write('File {} written\n'.format(f_name))

    # write list of uniques
    unique_occupancy = np.sum(consensus1[start_time:end_time, :], axis=0) / (end_time - start_time)
    name, ext = os.path.splitext(list_file)
    unique_list_file = name + '_unique' + ext
    backup_if_exists(unique_list_file)
    with open(unique_list_file, 'w') as f_handle:
        f_handle.write('trajectory: {}\n'.format(trajectory_file))
        f_handle.write('topology: {}\n'.format(topology_file))
        if args.dt:
            f_handle.write('dt: {}\n'.format(args.dt))
        f_handle.write('list_file: {}\n'.format(list_file))
        f_handle.write('group_file: {}\n'.format(groups_file))
        for pair, the_occupancy in zip(unique_pairs, unique_occupancy):
            label = pair[0].split('-')
            f_handle.write("{:>7} - {:<7} {:>3}%\n".format(label[0], label[1], int(round(the_occupancy * 100.0, 0))))
