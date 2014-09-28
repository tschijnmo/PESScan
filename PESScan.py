#!/usr/bin/env python

"""A simple script for generating configurations in a PES scan

The input file is going to be a simple JSON file specifying the configurations
that is needed for the program to run. The name of the JSON file should be given
as a command line argument. Specifically, it needs to have the following fields:

static-mol
  The file name for the coordinate of the molecule that is placed static in the
  PES scan. The format of the file should be a simple one of lines of element
  symbol followed by the three Cartesian coordinates. All blank lines and lines
  starting with the hash sign are going to be ignored.

scan-mol
  The file name for the coordinate of the molecule that is to be moved in the
  PES scan. It is going to be given in the same format as the static molecule.

translation-vector
  The translation vector to translate the molecule to scan.

from-point, to-pint
  Alternative to the translation vector. The molecule are going to be translated
  by the vector from the point given to the other point. If more than one points
  are given in a list, the centre of them are going to be used for the point.

translations
  A list of translations. Their product with the normalized translation vector
  are going to be used for translating the molecule.

even-mesh
  An alternative to the translation list, an even mesh of grid points are going
  to be generated. It is given as three numbers, with the first two being the
  initial and final distance to scan. The last number is going to be taken as
  the number of grid points between the initial and the final.

files
  A list of file names that is going to be configured by the code. They are
  going as python string template templates. The place holders that is supported
  are

  coords
    The coordinates of the molecules.

  dist
    The distance that is translated

  sn
    A one-based serial number for the data point

After the invocation of the program, a series of directories of integral names
are created for the data points. With all the files in the ``files`` list
configured and copied into each of them.

"""

import itertools
import sys
import json
import string
import math
import collections
import os
from os import path
import argparse

# pylint: disable=star-args


#
# The input reader
# ----------------
#

def get_input(file_name):

    """Gets the data structure from the JSON input file"""

    with open(file_name, 'r') as inp_file:
        inp = json.load(inp_file)
    return inp

def read_mol(file_name):

    """Reads a list of atoms from a file with file name

    The atoms are going to be returned as a quadruple of element symbol and
    Cartesian coordinates

    """

    with open(file_name, 'r') as in_f:

        raw_lines = [i.strip() for i in in_f]
        lines = [i for i in raw_lines if len(i) > 0 and i[0] != '#']

        atms = []
        for i in lines:
            fields = i.split()
            try:
                atms.append(
                    (fields[0], ) + tuple(float(j) for j in fields[1:4])
                    )
            except (IndexError, ValueError):
                print 'Corrupt input line in file %s:' % file_name
                print i
                raise IOError()
            continue

        return atms

def normalize(vec):

    """Normalizes a vector"""

    vec_list = list(vec)

    norm = math.sqrt(
        sum(i ** 2 for i in vec_list)
        )

    return tuple(i / norm for i in vec_list)

def centre(vecs):

    """Computes the centre of a list of points"""

    return tuple(
        sum(i) / len(vecs) for i in itertools.izip(*vecs)
        )

def get_transl_vec(inp):

    """Gets the translational vector from the input data structure"""

    if 'translation-vector' in inp:
        raw_vec = inp['translation-vector']
        return normalize(raw_vec)
    elif 'from-point' in inp and 'to-point' in inp:
        points = []
        for i in ['from-point', 'to-point']:
            raw_p = inp[i]
            if isinstance(raw_p[0], list):
                points.append(
                    centre(raw_p)
                    )
            else:
                points.append(raw_p)
        return normalize(
            j - i for i, j in itertools.izip(*points)
            )
    else:
        print "Translation vector is not given in the input!"
        raise IOError()

def get_grid(inp):

    """Gets the grid of data points to scan"""

    if 'translations' in inp:
        return inp['translations']
    elif 'even-mesh' in inp:
        mesh = inp['even-mesh']
        diff = (mesh[1] - mesh[0]) / (mesh[2] - 1)
        return [
            mesh[0] + diff * i for i in xrange(0, mesh[2] + 1)
            ]
    else:
        print 'Translation data point mesh is not given in the input!'
        raise IOError()


#
# Molecules translation
# ---------------------
#

def transl_mol(mol, vec):

    """Translates a molecule"""

    atms = []

    for atm_i in mol:
        atms.append(
            (atm_i[0], ) + tuple(i + j
                                 for i, j in itertools.izip(atm_i[1:4], vec))
            )
        continue

    return atms


#
# Scan data structure
# -------------------
#
# A scan point is going to be stored in a simple data structure having the atoms
# coordinates in field ``coords`` and the translated distance in ``dist``. The
# serial number is going to be read from its position in a list.
#

ScanPoint = collections.namedtuple(
    'ScanPoint',
    ['coords', 'dist']
    )

def gen_sp(static_mol, scan_mol, transl_vec, dist):

    """Generates a scan point data structure"""

    vec = tuple(i * dist for i in transl_vec)
    transled_scan_mol = transl_mol(scan_mol, vec)
    atms = static_mol + transled_scan_mol

    return ScanPoint(coords=atms, dist=dist)


#
# Dump a scan point to the file system
# ------------------------------------
#

def dump_sp(sp, file_names, dir_name):

    """Dumps a scan point into a directory

    :param sp: The scan point
    :param file_names: A list of file names to process
    :param dir_name: The directory name to dump the files

    """

    try:
        os.mkdir(dir_name)
    except OSError:
        pass

    atm_strs = [' %s   %f  %f  %f\n' % i for i in sp.coords]
    atm_strs[-1] = atm_strs[-1][0:-1]
    atm_str = ''.join(atm_strs)

    dist_str = '%f' % sp.dist
    sn_str = dir_name

    for file_name_i in file_names:
        with open(file_name_i, 'r') as file_i:

            content = file_i.read()
            template = string.Template(content)
            res = template.substitute(
                coords=atm_str,
                dist=dist_str,
                sn=sn_str
                )

            out_file_name = path.join(os.curdir, dir_name, file_name_i)
            with open(out_file_name, 'w') as out_file:
                out_file.write(res)

        continue

    return 0


#
# Main driver
# -----------
#

def main():

    """The main driver subroutine"""

    cur_module = sys.modules[__name__]
    parser = argparse.ArgumentParser(description=cur_module.__doc__)
    parser.add_argument('inp', help='The input file in JSON format')
    args = parser.parse_args()
    inp = get_input(args.inp)

    static_mol = read_mol(inp['static-mol'])
    scan_mol = read_mol(inp['scan-mol'])
    transl_vec = get_transl_vec(inp)
    grid = get_grid(inp)

    sps = [gen_sp(static_mol, scan_mol, transl_vec, i) for i in grid]

    dir_names_len = len(str(len(sps) + 1))
    dir_format = '%%%d.%dd' % (dir_names_len, dir_names_len)

    for i, v in enumerate(sps):
        dir_name = dir_format % (i + 1)
        dump_sp(v, inp['files'], dir_name)
        continue

    return 0


if __name__ == '__main__':
    main()
