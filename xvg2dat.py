#!/usr/bin/env python

from __future__ import print_function, division
import sys
import numpy as np
import os
import argparse
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


def main():
    ######################
    #                    #
    # Parse command line #
    #                    #
    ######################
    p = argparse.ArgumentParser()
    p.add_argument('-dir', required=False, dest='directory', action=check_write_dir(),
                   help='work directory')
    args = p.parse_args()
    directory = os.getcwd()
    if args.directory:
        directory = os.path.abspath(args.directory)
    files = os.listdir(directory)
    for file in files:
        name, ext = os.path.splitext(file)
        if ext != '.xvg':
            continue
        datname = os.path.abspath(name) + '.dat'
        with open(file, 'r') as handle1:
            _, data = read_xvg(handle1)
            backup_if_exists(datname)
            with open(datname,"w") as handle2:
                for i, datum in enumerate(data[1]):
                    handle2.write('{} {}\n'.format(i, float(datum)*10))


if __name__ == '__main__':
    main()
