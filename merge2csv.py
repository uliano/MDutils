#!/usr/bin/env python

# old file its function should now be part of do_hbonds

from __future__ import print_function, division
import sys
import numpy as np
import os
import argparse

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


def main():
    ######################
    #                    #
    # Parse command line #
    #                    #
    ######################
    p = argparse.ArgumentParser()
    p.add_argument('-f1', required=True, dest='file1', action=check_write_ext('csv'),
                   help='first .csv will be replaced with the merged')
    p.add_argument('-f2', required=True, dest='file2', action=check_write_ext('csv'),
                   help='second .csv')
    args = p.parse_args()
    with open(args.file1, 'r') as file:
        lines1 = file.readlines()
    with open(args.file2, 'r') as file:
        lines2 = file.readlines()
    os.rename(args.file1, args.file1+'.merged')
    os.rename(args.file2, args.file2+'.merged')
    file = open(args.file1, "w")
    header = lines1.pop(0)
    file.write(header)
    lines2.pop(0)
    for i in range(len(lines1)):
        time, dist1, occ1 = lines1[i].strip().split(',')
        _, dist2, occ2 = lines2[i].strip().split(',')
        dist = str(min(float(dist1), float(dist2)))
        occ = 'None'
        if occ1 == 'Present' or occ2 == 'Present':
            occ = 'Present'
        file.write('{},{},{}\n'.format(time, dist, occ))
    file.close()




if __name__ == '__main__':
    main()
