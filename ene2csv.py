#!/usr/bin/env python

# from a bunch of alphabetically ordered CHARMM .ENE files creates a single .csv

from __future__ import print_function, division
from functools import reduce
import sys
import os
import argparse


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


def read_header(lines):
    # lots of heuristics here as no file format specification was available
    first = 1
    last = first
    while lines[last][1].isupper():
        last += 1
    nlines = last - first
    header_labels = ""
    for line in lines[first:first+nlines]:
        header_labels += line
    column_names = header_labels.split()
    fields = [reduce((lambda x, y: x + y), lines[i: i + nlines]).split() for i in range(1 + nlines, len(lines), nlines)]
    return column_names, fields




def main():
    ######################
    #                    #
    # Parse command line #
    #                    #
    ######################
    p = argparse.ArgumentParser()
    p.add_argument('-dir', required=False, dest='directory', action=check_write_dir(),
                   help='work directory', default=os.getcwd())
    p.add_argument('-o', required=False, dest='output', action=check_write_dir(),
                   help='work directory', default='energy.csv')
    args = p.parse_args()
    directory = os.getcwd()
    if args.directory:
        directory = os.path.abspath(args.directory)
    files = os.listdir(directory)
    datname = os.path.join(directory, args.output)
    write_header = True
    backup_if_exists(datname)
    with open(datname, "w") as handle2:
        for file in sorted(files):
            name, ext = os.path.splitext(os.path.basename(file))
            if ext != '.ENE':
                continue
            print("reading " + file)
            with open(file, 'r') as handle1:
                lines = handle1.readlines()
                lines = [line.rstrip() for line in lines]
                header, rows = read_header(lines)
                if write_header:
                    handle2.write(",".join(header)+"\n")
                    write_header = False
                for row in rows:
                    handle2.write(",".join(row)+"\n")


if __name__ == '__main__':
    main()
