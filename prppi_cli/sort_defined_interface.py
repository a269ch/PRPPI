#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Module of general"""

import json
import os
import click
import time
import warnings
from optparse import OptionParser

from Bio.PDB import PDBExceptions, PDBParser
from Bio.PDB.Model import Model
from numpy import float32

__author__ = 'Aleksei Che'
__version__ = '0.1.0'

from .settings import PATH_TO_RESULT, standard_aa_names


def min_distance(pose: Model, residue1: tuple, residue2: type) -> float32:
    """
    :param pose:
    :param residue1:
    :param residue2:
    :return: The minimal distance between any pair of atoms between two residues
    """
    min_dist = -1
    chain1, seqpos1, resid1 = residue1
    chain2, seqpos2, resid2 = residue2
    for atom1 in pose[chain1][int(seqpos1)]:
        if 'H' in atom1.get_name():
            continue
        for atom2 in pose[chain2][int(seqpos2)]:
            if 'H' in atom2.get_name():
                continue
            if atom1 - atom2 < min_dist or min_dist == -1:
                min_dist = atom1 - atom2
    return min_dist


def get_residue_sides(pose: Model, side1_chains: str, side2_chains: str) -> tuple:
    """
    :param pose:
    :param side1_chains:
    :param side2_chains:
    :return: Two lists with residues in a complex split by their interface side
    """
    side1_residues = []
    for chain in side1_chains:
        for res in pose[chain]:
            seqpos = res.get_id()[1]
            if res.get_resname() in standard_aa_names:
                resid = standard_aa_names[res.get_resname()]
                side1_residues.append((chain, str(seqpos), resid))

    side2_residues = []
    for chain in side2_chains:
        for res in pose[chain]:
            seqpos = res.get_id()[1]
            if res.get_resname() in standard_aa_names:
                resid = standard_aa_names[res.get_resname()]
                side2_residues.append((chain, str(seqpos), resid))
    return side1_residues, side2_residues


def sc_distance_cutoff(
        pose: Model,
        side1_chains: str,
        side2_chains: str,
        cutoff: float,
        groups: dict
) -> dict:
    """
    :param pose:
    :param side1_chains:
    :param side2_chains:
    :param cutoff
    :param groups
    :return: List of residues with an atom within 5 A of an atom on the opposing chain
    """
    sides = {}
    side1_residues, side2_residues = get_residue_sides(pose, side1_chains, side2_chains)

    for res1 in side1_residues:
        key_amino = ''.join(res1)
        nearby_atom_cutoff = [k for k, v in groups.items() if key_amino[-1] in v][0] if groups else cutoff

        sides_dict = {
            ''.join(res1_self): distance for res1_self in side1_residues
            if res1 != res1_self and (distance := float(min_distance(pose, res1, res1_self))) <= nearby_atom_cutoff
        }

        sides_dict2 = {
            ''.join(res1_self): distance for res1_self in side2_residues
            if res1 != res1_self and (distance := float(min_distance(pose, res1, res1_self))) <= nearby_atom_cutoff
        }
        sides[key_amino] = {side1_chains: sides_dict, side2_chains: sides_dict2}

        if not sides[key_amino][side2_chains]:
            sides.pop(key_amino)

        else:
            res2_b = min([v for k, v in sides_dict2.items() if side2_chains in k])  # minimum
            res2_b_1 = [key for key in sides_dict2 if sides_dict2[key] == res2_b][0]  # key minimum side 2
            res2 = res2_b_1[0], res2_b_1[1:-1], res2_b_1[-1]

            sides2_dict1 = {
                ''.join(res1_self): distance_s2 for res1_self in side1_residues if res1 != res1_self
                and (distance_s2 := float(min_distance(pose, res2, res1_self))) <= cutoff
            }

            sides2_dict2 = {
                ''.join(res1_self): distance_s2 for res1_self in side2_residues if res1 != res1_self
                and (distance_s2 := float(min_distance(pose, res2, res1_self))) <= cutoff
            }

            sides[key_amino][side1_chains] = {
                amino[:-1]: sides_dict[amino] for amino in sides_dict if amino in sides2_dict1
            }

            sides[key_amino][side2_chains] = {
                amino[:-1]: sides_dict2[amino] for amino in sides_dict2 if amino in sides2_dict2
            }

    return sides


def run(pdb, side_1, side_2, cutoff, groups, out_file='out'):
    """
    :return:
    """
    parser = PDBParser()
    structure = parser.get_structure('X', pdb)
    # model: Model = structure[0]
    # print(model['A'].get_unpacked_list()[0])
    time.sleep(10)
    result = sc_distance_cutoff(structure[0], side_1, side_2, cutoff, groups)
    path = f"{PATH_TO_RESULT}"
    try:
        os.mkdir(path)
    except OSError:
        pass

    with open(f"{path}/{out_file}.json", "w") as write_file:
        json.dump(result, write_file, indent=4)
        return write_file.name


if __name__ == '__main__':
    usage = "%prog [options] <pdb_file>"
    parser = OptionParser(usage)
    parser.add_option("--output", dest="output", default="out", help="Output name for json file")
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error('specify a pdb file or use -h')
    elif len(args) > 1:
        print('Warning: only the first pdb is considered by this script. The rest will be ignored')

    warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)

    print('Processing', args[0])
