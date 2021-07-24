import json
import pprint

import click

from . import sort_defined_interface
from .osprey_algorithms import bbkstar
from .settings import SIDE_1, SIDE_2

GROUPS = None


@click.command()
@click.argument('pdb', nargs=1, required=True)
@click.option('--side_1', '-s1', default=SIDE_1, show_default=True)
@click.option('--side_2', '-s2', default=SIDE_2, show_default=True)
@click.option('--cutoff', default=5.0, show_default=True)
@click.option('--groups', default=GROUPS, show_default=True)
# @click.option('--algorithm', type=click.Choice(['bbkstar', 'kstar', 'markstar_bbkstar', 'markstar_kstar']))
def cli(pdb, side_1, side_2, cutoff, groups):
    """ Program for protein redesign to increase affinity of protein protein interaction. """
    sorted_json = sort_defined_interface.run(pdb, side_1, side_2, cutoff, groups)
    input('The calculation of the data is finished, press Enter to see the result')
    click.edit(filename=sorted_json, editor='nano')
    input('press Enter to run OSPREY')

    with open(sorted_json) as sj:
        sorted_amino = json.load(sj)

    with click.progressbar(sorted_amino.keys()) as bar:
        for amino in bar:
            bbkstar.run(
                pdb=pdb,
                interest_aa=amino[:-1],
                neighbours_side1=list(sorted_amino.get(amino).get(side_1).keys()),
                neighbours_side2=list(sorted_amino.get(amino).get(side_2).keys()),
            )


if __name__ == '__main__':
    cli()
