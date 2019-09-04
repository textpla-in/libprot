# -*- coding: utf-8 -*-
import re
import sys
import typing

import click
import numpy as np
import requests
from Bio import SeqIO
from Bio.PDB import PDBParser, Superimposer, PPBuilder
from Bio.SeqRecord import SeqRecord

from libprot.pdb import get_amino_acids
from libprot import amber
from libprot.amber import ForceFieldType, reduce_pdb


def error(msg):
    print(msg, file=sys.stderr)
    sys.exit(1)


@click.group()
def main():
    pass


@click.command()
@click.option('-i', '--id', 'pdb_id', required=True)
def fetch_pdb(pdb_id):
    pdb_file = f'{pdb_id}.pdb'
    url = f'https://files.rcsb.org/download/{pdb_file}'
    r = requests.get(url)
    with open(pdb_file, 'w') as f:
        f.write(r.text)

    print('saved', pdb_file)


@click.command()
@click.argument('pdb1', type=click.Path(exists=True))
@click.argument('pdb2', type=click.Path(exists=True))
def calc_rmsd(pdb1, pdb2):
    parser = PDBParser()
    first = parser.get_structure(pdb1.replace('.pdb', ''), pdb1)
    second = parser.get_structure(pdb2.replace('.pdb', ''), pdb2)

    fixed = list(first.get_atoms())
    moving = list(second.get_atoms())
    if len(fixed) != len(moving):
        error(f'There is a different number of atoms in {pdb1} ({len(fixed)}) than in {pdb2} ({len(moving)})')

    sup = Superimposer()
    # fixed and moving are lists of atom objects
    sup.set_atoms(fixed, moving)


def edit_distance(seq1, seq2):
    mat = np.zeros(shape=(len(seq1) + 1, len(seq2) + 1), dtype=int)

    # base cases
    for idx in range(mat.shape[0]):
        mat[idx, 0] = idx

    for idx in range(mat.shape[1]):
        mat[0, idx] = idx

    def transition_fn(row, col):
        from_diag = mat[row - 1, col - 1] + int(seq1[row - 1] != seq2[col - 1])
        from_above = mat[row - 1, col] + 1
        from_left = mat[row, col - 1] + 1

        mat[row, col] = min(from_diag, from_above, from_left)

    indices = [(row, col) for row in range(1, mat.shape[0]) for col in range(1, mat.shape[1])]
    for row, col in indices:
        transition_fn(row, col)

    return mat[len(seq1), len(seq2)]


def to_fasta_impl(structure, pdb_id, description, mutations, out_stream):
    if pdb_id is None:
        pdb_id = ""
    if description is None:
        description = ""

    ppb = PPBuilder()
    seq = [pp.get_sequence() for pp in ppb.build_peptides(structure)][0]

    if mutations:
        mut_seq = seq.tomutable()
        mutations = [re.search(r'(?P<original>\D)(?P<index>\d+)(?P<mutation>\D)', line) for line in mutations]
        for mut in mutations:
            idx = int(mut.group('index')) - 1
            old_aa = mut.group('original')
            new_aa = mut.group('mutation')
            if old_aa != mut_seq[idx]:
                raise Exception(
                    f'The amino acid in {pdb_id} at index {mut.group("index")} is {mut_seq[idx]}, not {old_aa}')
            mut_seq[idx] = new_aa
        seq = mut_seq.toseq()

    sequence_record = SeqRecord(seq, id=pdb_id, description=description)
    SeqIO.write([sequence_record], out_stream, "fasta")


@click.command()
@click.argument('pdb', type=click.Path(exists=True))
@click.option('-i', '--id', 'pdb_id', required=False)
@click.option('-d', '--description', 'description', required=False)
@click.option('-m', '--mutation-file', 'mut_file', type=click.Path(exists=True), required=False)
def to_fasta(pdb, pdb_id, description, mut_file):
    mutations = []
    if mut_file:
        with open(mut_file, 'r') as mut:
            mutations = [line.strip() for line in mut.readlines()]

    parser = PDBParser()
    structure = parser.get_structure(pdb.replace('.pdb', ''), pdb)

    to_fasta_impl(structure, pdb_id, description, mutations, sys.stdout)


@click.command()
@click.argument('path', type=click.Path(exists=True))
def dump_residues(path):
    amino_acids = get_amino_acids(path)
    for aa in amino_acids:
        print(aa)


def file_or_stdin(f):
    return f or '-'


def print_command_output(return_code, stdout, stderr):
    if return_code:
        print('Reduce returned a non-zero error code.', file=sys.stderr)
        print(stderr.read(), file=sys.stderr)
    else:
        print(stdout.read())


@click.command(help="Run reduce on the PDB file or stdin")
@click.option('--file', '-f', type=click.Path(exists=True))
def run_reduce_on_pdb(file):
    with click.open_file(file_or_stdin(file)) as f:
        print_command_output(*reduce_pdb(f))


@click.command(help="Run pdb4amber on the PDB file or stdin")
@click.option('--file', '-f', type=click.Path(exists=True, allow_dash=True))
def run_pdb4amber(file):
    with click.open_file(file_or_stdin(file)) as f:
        print_command_output(*amber.prep_pdb_for_amber(f))


__ff_map = {
    'protein': ForceFieldType.PROTEIN,
    'dna': ForceFieldType.DNA,
    'rna': ForceFieldType.RNA,
    'water': ForceFieldType.WATER,
    'general': ForceFieldType.GENERAL_FF
}


@click.command(help="Use amber to parameterize a PDB file or stdin")
@click.option('--file', '-f', type=click.Path(exists=True, allow_dash=True))
@click.option('--molecule-type', type=click.Choice(__ff_map.keys()))
def parameterize_pdb_file(file, molecule_type):
    forcefield_type = __ff_map[molecule_type]
    with click.open_file(file_or_stdin(file)) as f:
        return_code, prmtop, prmcrd = amber.make_parameter_files(f, forcefield_type)
        print(f'prmtop\n========')
        print(prmtop.read())
        print(f'prmcrd\n========')
        print(prmcrd.read())


@click.command(help="Use sander to minimize a molecule")
@click.option('--topology', '-t', type=click.File())
@click.option('--coordinates', '-c', type=click.File())
@click.option('--rounds', '-c', type=click.IntRange(min=100, max=1000))
def minimize_parameter_files(topology, coordinates, rounds):
    print_command_output(*amber.run_sander(topology, coordinates, rounds))


@click.command(help="Convert an amber coordinates file to a PDB")
@click.option('--topology', '-t', type=click.File(mode='r'))
@click.option('--coordinates', '-c', type=click.File(mode='rb'))
def coordinates_to_pdb(topology: typing.TextIO, coordinates: typing.BinaryIO):
    print_command_output(*amber.convert_coords_to_pdb(topology, coordinates))

@click.command(help="Minimize a protein in a PDB file")
@click.option('--file', '-f', type=click.File())
def minimize_protein(file):
    print_command_output(*amber.do_minimize_pdb(file, ForceFieldType.PROTEIN, 200))


main.add_command(fetch_pdb)
main.add_command(calc_rmsd)
main.add_command(to_fasta)
main.add_command(dump_residues)
main.add_command(run_reduce_on_pdb, name='run-reduce')
main.add_command(run_pdb4amber)
main.add_command(parameterize_pdb_file, name='parameterize')
main.add_command(coordinates_to_pdb, name='coords-to-pdb')
main.add_command(minimize_parameter_files, name='run-sander')
main.add_command(minimize_protein, name='minimize-protein')

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
