#!/usr/bin/env python

# probably an older version of do_hbonds

from __future__ import print_function, division
import sys
import numpy as np
import pandas as pd
import os
import argparse


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


def process_meta(meta):
    uncommented = []
    for i in meta:
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
        meta = [fhandle.readline()]
        while not meta[-1].startswith("static char *gromacs_xpm[]"):
            meta.append(fhandle.readline())

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
                meta.append(i)
                continue
            j = unquote(i)
            colors.append([lut[j[k:k + nb]][0] for k in range(0, nx, nb)])
            values.append([lut[j[k:k + nb]][1] for k in range(0, nx, nb)])
        meta = process_meta(meta)

        if reverse:
            values.reverse()
            colors.reverse()

        return colors, values, meta, lut


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
    lines = fhandle.readlines()
    lines.pop(0)  # remove comment
    n_atoms = int(lines.pop(0))
    lines.pop(-1)  # remove last line containing the box
    # box = [float(i) for i in box_list]
    atoms = dict()
    residues = dict()
    res_adder = 0
    atom_adder = 0
    res_step = 100000
    atom_step = 100000
    for line in lines:
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
        atoms[a_num] = {'aname': a_nam, 'rname': r_nam, 'rnum': r_num}
        if r_num in residues:
            residues[r_num]['atoms'].append(a_num)
        else:
            residues[r_num] = {'rname': r_nam, 'atoms': [a_num]}
    assert n_atoms == len(atoms)
    return atoms, residues


def read_pdb(fhandle):
    atoms = dict()
    residues = dict()
    res_adder = 0
    atom_adder = 0
    res_step = 10000
    atom_step = 100000

    for line in fhandle:
        if line[:4] == 'ATOM' or line[:6] == "HETATM":
            r_nam = line[17:21].strip()
            rn = line[22:26]
            r_num = str(int(rn) + res_adder)
            a_nam = line[12:16].strip()
            an = line[6:11]
            a_num = str(int(an) + atom_adder)
            if rn == '9999':
                res_adder += res_step
            if an == '99999':
                atom_adder += atom_step
            atoms[a_num] = {'aname': a_nam, 'rname': r_nam, 'rnum': r_num}
            if r_num in residues:
                residues[r_num]['atoms'].append(a_num)
            else:
                residues[r_num] = {'rname': r_nam, 'atoms': [a_num]}
    return atoms, residues


def main():
    indexes = None

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=argparse.FileType('r'), required=True,
                        help='gromacs hbond index file .ndx', dest='index_file')
    parser.add_argument('-s', type=argparse.FileType('r'), required=True,
                        help='topology file .gro', dest='topology_file')
    parser.add_argument('-m', type=argparse.FileType('r'), required=True,
                        help='gromacs hbond matrix file .xpm', dest='matrix_file')
    parser.add_argument('-l', nargs='?', type=argparse.FileType('w'), required=False,
                        help='list output file', dest='list_file', default='hbond_list.txt')
    parser.add_argument('-o', nargs='?', type=argparse.FileType('w'), required=False,
                        help='index file for gmx distance, defaults to distance.ndx', dest='output_file',
                        default='distance.ndx')
    parser.add_argument('-c', nargs='?', type=argparse.FileType('w'), required=False,
                        help='hbond consensus, defaults to consensus.csv', dest='consensus_file',
                        default='consensus.csv')
    parser.add_argument('-b', nargs='?', type=int, dest='begin', required=False, help='begin time (ps)', default=0)
    parser.add_argument('-e', nargs='?', dest='end', required=False, help='end time (ps)', default='No')
    parser.add_argument('-t', nargs='?', type=int, dest='threshold', required=False,
                        help='threshold presence percentage to accept hbond, defaults to 0',
                        default=0)
    parser.add_argument('-sel', nargs='?', dest='select_file', required=False, type=argparse.FileType("w"),
                        help='Selection file for gmx distance, defaults to same name of the OUTPUT_FILE'
                             ' but with .sel extension')
    parser.add_argument('-min-res-dist', nargs='?', dest='min_res_dist', required=False, type=int, default=0,
                        help='Consider hbonds with residue number difference greater than, defults to 0 ')
    args = parser.parse_args()
    parser.parse_args
    matrix_name = args.matrix_file.name
    index_name = args.index_file.name
    topology_name = args.topology_file.name

    if args.index_file.name[-4:] == '.ndx':
        indexes = read_index(args.index_file)
    else:
        sys.stderr.write("-f requires a .ndx file.\n")
        exit(2)
    args.index_file.close()
    if args.topology_file.name[-4:] == '.gro':
        atoms, _ = read_gro(args.topology_file)
    elif args.topology_file.name[-4:] == '.pdb':
        atoms, _ = read_pdb(args.topology_file)
    else:
        sys.stderr.write("-s requires a .pdb or a .gro file.\n")
        sys.exit(2)
    args.topology_file.close()
    if args.matrix_file.name[-4:] == '.xpm':
        _, matrix, meta, _ = read_xpm(args.matrix_file, reverse=True)
    else:
        sys.stderr.write("-m requires a .xpm file.\n")
        sys.exit(2)
    args.matrix_file.close()
    if not args.output_file:
        sys.stderr.write("Unable to open output file.")

    threshold = float(args.threshold) / 100.0

    x_values = meta['x-axis'].split()
    x_values = [int(i) for i in x_values]

    if args.end == 'No':
        end = x_values[-1]
    else:
        end = args.end

    first_column = 0
    last_column = len(x_values) - 1
    for i in range(len(x_values)):
        if x_values[i] >= args.begin:
            first_column = i
            break
    for i in reversed(range(len(x_values))):
        if x_values[i] <= end:
            last_column = i
            break

    _, lines = indexes[-1]

    matrix = np.array(matrix)
    vector = np.sum(matrix[:, first_column:last_column], axis=1) / (last_column - first_column)

    args.list_file.write("Begin = {0}\nend = {1}\nthreshold = {2} %\n".format(x_values[first_column],
                                                                              x_values[last_column], args.threshold))
    args.list_file.write("Index file: {0}\nmatrix file: {1}\ntopology: {2}\n".format(index_name,
                                                                                     matrix_name,
                                                                                     topology_name))
    args.list_file.write("Minimal difference between residues = {}\n\n".format(args.min_res_dist))

    args.list_file.write("  #  |           Donor         |   Hydrogen  |         Acceptor        | %   \n")
    args.list_file.write("---------------------------------------------------------------------------\n")

    if args.select_file:
        select_file = args.select_file
    else:
        select_file = open(str(os.path.abspath(args.output_file.name))[:-4]+".sel", "w")
    df = pd.DataFrame()
    j = 0
    i = 0
    for line, value in zip(lines, vector):
        if value > threshold:
            dnum, hnum, anum = line.strip().split()
            dname = atoms[dnum]['aname']
            drname = atoms[dnum]['rname']
            drnum = atoms[dnum]['rnum']
            hname = atoms[hnum]['aname']
            aname = atoms[anum]['aname']
            arname = atoms[anum]['rname']
            arnum = atoms[anum]['rnum']
            if abs(int(drnum) - int(arnum)) > args.min_res_dist:
                args.list_file.write("{:>4} | ".format(str(i)))
                args.list_file.write("{:>4} {:6} {:>4} {:6} | ".format(drname, drnum, dname, dnum))
                args.list_file.write("{:>4} {:6} | ".format(hname, hnum))
                args.list_file.write("{:>4} {:6} {:>4} {:6} | {:>2}\n".format(arname, arnum, aname, anum,
                                                                              int(round(value * 100.0, 0))))
                # args.output_file.write("[ {}-{}_{}-{}<->{}-{}_{}-{}_({}_{}) ]\n".format(drname, drnum, dname,
                # dnum, arname, arnum, aname,
                # anum, j, i))
                df["{}{}_{}{}_{}".format(drname, drnum, arname, arnum, i)] = matrix[j, :]
                args.output_file.write("[ {}{}_{}{}_{} ]\n".format(drname, drnum, arname, arnum, i))
                args.output_file.write("{0:>6} {1:>6}\n".format(dnum, anum))
                select_file.write(str(i))
                select_file.write("\n")
                i += 1
        j += 1
    df.to_csv(args.consensus_file)


    args.list_file.close()
    args.output_file.close()
    select_file.close()


if __name__ == '__main__':
    main()
