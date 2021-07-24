#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from os import path
path.abspath(__file__)
from yaml import load

try:
    from yaml import CDumper as Dumper
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Dumper, Loader

CURRENT_DIR = path.dirname(path.abspath(__file__))

with open(f'{CURRENT_DIR}/settings.yml', 'r') as f:
    config = load(f, Loader=Loader)

standard_aa_names = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K",
    "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y", "TYS": "Y"
}


PATH: dict = config.get('PATH')
PATH_TO_RESULT = f"{PATH.get('PATH_TO_RESULT')}result"

CUTOFFS: dict = config.get('CUTOFFS')
GROUPS = CUTOFFS.get('GROUPS')
GROUPS_SIDE1 = GROUPS.get('SIDE1')
GROUPS_SIDE2 = GROUPS.get('SIDE2')

SIDES: dict = config.get('SIDES')
SIDE_1 = SIDES.get('SIDE_1')
SIDE_2 = SIDES.get('SIDE_2')

RESIDUES_SIDE1: dict = config.get('RESIDUES_SIDE1')
RESIDUES_SIDE2: dict = config.get('RESIDUES_SIDE2')

OSPREY: dict = config.get('OSPREY')
CPU = OSPREY.get('CPU')
GPU = OSPREY.get('GPU')
RAM = OSPREY.get('RAM')
EPSILON = OSPREY.get('EPSILON')
MUTATIONS = OSPREY.get('MUTATIONS')
