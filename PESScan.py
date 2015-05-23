#!/usr/bin/env python3

"""A simple script for generating configurations in a PES scan

The input file is going to be a simple YAML file specifying the
configurations that is needed for the program to run. The name of the YAML
file should be given as a command line argument. Specifically, it needs to
have the following fields:

static-mol
  The file name for the coordinate of the molecule that is placed static in
  the PES scan. The format of the file should be a simple one of lines of
  element symbol followed by the three Cartesian coordinates. All blank lines
  and lines starting with the hash sign are going to be ignored.

scan-mol
  The file name for the coordinate of the molecule that is to be moved in the
  PES scan. It is going to be given in the same format as the static molecule.

translation-vector
  The translation vector to translate the molecule to scan.

from-point, to-point
  Alternative to the translation vector. The molecule are going to be
  translated by the vector from the point given to the other point. If more
  than one points are given in a list, the centre of them are going to be
  used for the point.

from-atom, to-atom
  Yet another alternative for giving the translation vector. Here the
  coordinates are given as indices of the atoms in the molecules. The
  ``from-atom`` is based on the static molecules, while to ``to-atom`` is
  based on the scan molecule. It can be both a single integer or a list of
  integers. When a list of given, the centre of the atoms are going to be
  used. The indices are all one-based.

translations
  A list of translations. Their product with the normalized translation
  vector are going to be used for translating the molecule.

even-mesh
  An alternative to the translation list, an even mesh of grid points are
  going to be generated. It is given as three numbers, with the first two
  being the initial and final distance to scan. The last number is going to
  be taken as the number of grid points between the initial and the final.

files

  A list of file names that is going to be configured by the code. They are
  going to be names of Jinja templates in the current working directory. The
  fields that are supported are

  atoms

      The atoms in the molecules, given as a dictionary with keys ``element``
      for the element and ``coords`` for the Cartesian coordinates. Also
      ``idx`` is given for the index of the atom in the list of all atoms.

  translation
    The distance that is translated

  distance
    The distance between the from and to points after the translation.

  sn
    A one-based serial number for the data point

After the invocation of the program, a series of directories of integral names
are created for the data points. With all the files in the ``files`` list
configured and copied into each of them.

"""


import sys
import math
import collections
import os
from os import path
import argparse

import yaml
import jinja2


# pylint: disable=star-args


#
# The input reader
# ----------------
#


def get_input(file_name):
    """Gets the data structure from the YAML input file"""

    with open(file_name, 'r') as inp_file:
        inp = yaml.load(inp_file)
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
                print('Corrupt input line in file {}:'.format(file_name))
                print(i)
                raise IOError()
            continue

        return atms


def normalize(vec, return_norm=False):
    """Normalizes a vector"""

    vec_list = list(vec)

    norm = math.sqrt(
        sum(i ** 2 for i in vec_list)
        )

    normed_vec = tuple(i / norm for i in vec_list)

    return (normed_vec, norm) if return_norm else normed_vec


def centre(vecs):
    """Computes the centre of a list of points"""

    return tuple(
        sum(i) / len(vecs) for i in zip(*vecs)
        )


def get_transl_vec(inp, static_mol, scan_mol):
    """Gets the translational vector from the input data structure

    Both the normalized translational vector and the original norm of the
    vector are going to be returned.
    """

    # pre-process the input object to support the indices of the atoms in the
    # molecules

    def idx2coord(idxes, mol):
        """Convert the indices to coordinates"""
        if isinstance(idxes, list):
            return [
                mol[i - 1][1:4] for i in idxes
                ]
        else:
            return mol[idxes - 1][1:4]
    if 'from-atom' in inp:
        inp['from-point'] = idx2coord(inp['from-atom'], static_mol)
    if 'to-atom' in inp:
        inp['to-point'] = idx2coord(inp['to-atom'], scan_mol)

    if 'translation-vector' in inp:
        raw_vec = inp['translation-vector']
        return normalize(raw_vec)
    elif 'from-point' in inp and 'to-point' in inp:
        points = []
        for i in ['from-point', 'to-point']:
            raw_p = inp[i]
            if isinstance(raw_p[0], list) or isinstance(raw_p[0], tuple):
                points.append(
                    centre(raw_p)
                    )
            else:
                points.append(raw_p)
        return normalize(
            (j - i for i, j in zip(*points)),
            return_norm=True
        )
    else:
        print("Translation vector is not given in the input!")
        raise IOError()


def get_grid(inp):
    """Gets the grid of data points to scan"""

    if 'translations' in inp:
        return inp['translations']
    elif 'even-mesh' in inp:
        mesh = inp['even-mesh']
        diff = (mesh[1] - mesh[0]) / (mesh[2] - 1)
        return [
            mesh[0] + diff * i for i in range(0, mesh[2] + 1)
            ]
    else:
        print('Translation data point mesh is not given in the input!')
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
            (atm_i[0], ) + tuple(i + j for i, j in zip(atm_i[1:4], vec))
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

def dump_sp(sp, orig_dist, file_names, dir_name): # pylint: disable=invalid-name
    """Dumps a scan point into a directory

    :param sp: The scan point
    :param orig_dist: The original distance between the from and to points.
    :param file_names: A list of file names to process
    :param dir_name: The directory name to dump the files
    """

    try:
        os.mkdir(dir_name)
    except OSError:
        pass

    atoms = [
        {
            'idx': i,
            'element': v[0]
            'coords': v[1:3]
        }
        for i, v in enumerate(sp.coords)
    ]

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(os.getcwd())
    )

    for file_name_i in file_names:
        templ = env.get_template(file_name_i)
        ctx = {
            'atoms': atoms,
            'translation': sp.dist,
            'distance': orig_dist + sp.dist,
            'sn': dir_name,
        }
        out_file_name = path.join(os.curdir, dir_name, file_name_i)
        templ.stream(ctx).dump(out_file_name)
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
    transl_vec, orig_dist = get_transl_vec(inp, static_mol, scan_mol)
    grid = get_grid(inp)

    sps = [gen_sp(static_mol, scan_mol, transl_vec, i) for i in grid]

    dir_names_len = len(str(len(sps) + 1))
    dir_format = '%%%d.%dd' % (dir_names_len, dir_names_len)

    for i, v in enumerate(sps): # pylint: disable=invalid-name
        dir_name = dir_format % (i + 1)
        dump_sp(v, orig_dist, inp['files'], dir_name)
        continue

    return 0


if __name__ == '__main__':
    main()
