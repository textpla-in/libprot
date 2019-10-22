import dataclasses
import itertools
import re
import sys
import typing
from pathlib import Path

import click
import numpy as np
import requests
from Bio.PDB import PDBParser, Superimposer, PPBuilder, Structure, is_aa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from libprot import amber
from libprot.amber import ForceFieldType, reduce_pdb
from libprot.pdb import get_amino_acids
from libprot.types import Mutation, Indel


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


def parse_mutations(mutations: typing.List[str]) -> typing.List[Mutation]:
    mutations = [re.search(r'(?P<chain>\D):(?P<mut_from>\D)(?P<index>\d+)(?P<mut_to>\D)', mutation) for mutation in mutations]
    return [Mutation(mut.group('chain'), int(mut.group('index')), mut.group('mut_from'), mut.group('mut_to')) for mut in mutations]


def index_in_lists(seqs: typing.List[typing.Any], idx: int) -> typing.Tuple[int, int]:
    seq_idx = 0
    for seq in seqs:
        try:
            seq[idx]
        except IndexError:
            idx -= len(seq)
            seq_idx += 1

    return seq_idx, idx


def mutate_residues(structure: Structure, seqs: typing.List[Seq], mutations: typing.List[Mutation]) -> typing.List[Seq]:

    def match_residue(res):
        return int(res.get_id()[1]) == mutation.index and res.get_parent().get_id() == mutation.chain

    for mutation in mutations:
        res = next(filter(match_residue, structure.get_residues()))
        index = [res for res in structure.get_residues() if is_aa(res)].index(res)
        seq_idx, idx_in_seq = index_in_lists(seqs, index)
        mut_seq = seqs[seq_idx].tomutable()

        if mutation.from_aa != mut_seq[idx_in_seq]:
            raise Exception(f'The amino acid at index {idx_in_seq} is {mut_seq[idx_in_seq]}, not {mutation.from_aa}')
        mut_seq[idx_in_seq] = mutation.to_aa
        seqs[seq_idx] = mut_seq.toseq()

    return seqs


def parse_to_structure(path: Path) -> Structure:
    return PDBParser(QUIET=True).get_structure(path.name, path)


def parse_indels(indels: typing.List[str]) -> typing.List[Indel]:
    indels = [re.search(r'(?P<chain>\D):(?P<amino_acid>\D)(?P<index>\d+)', insertion) for insertion in indels]
    indels = [Indel(ins.group('chain'), int(ins.group('index')), ins.group('amino_acid')) for ins in indels]
    return indels


def add_indels(structure: Structure, seqs: typing.List[Seq], indels: typing.List[Indel]) -> typing.List[Seq]:

    def match_residue(res):
        return int(res.get_id()[1]) == insertion.index and res.get_parent().get_id() == insertion.chain

    insertions = sorted(id for id in indels if id.is_insertion())
    deletions = sorted(id for id in indels if not id.is_insertion())

    while insertions:
        insertion, insertions = insertions[0], insertions[1:]
        res = next(filter(match_residue, structure.get_residues()))
        aas = [aa for aa in structure.get_residues() if is_aa(aa)]
        index = aas.index(res)
        seq_idx, idx_in_seq = index_in_lists(seqs, index)
        seq = seqs[seq_idx].tomutable()
        beginning, end = seq[:idx_in_seq], seq[idx_in_seq:]
        seqs[seq_idx] = (beginning + insertion.aa + end).toseq()

    return seqs


def extract_sequence_from_pdb_file(path: Path) -> typing.List[Seq]:
    return [pp.get_sequence() for pp in PPBuilder().build_peptides(parse_to_structure(path))]


def mutate_and_indels(structure: Structure, seqs: typing.List[Seq], mutations: typing.List[Mutation], indels: typing.List[Indel]) -> typing.List[Seq]:
    seqs = mutate_residues(structure, seqs, mutations)
    seqs = add_indels(structure, seqs, indels)
    return seqs


def to_fasta_impl(path: Path, pdb_id: str, description: str, mutations: typing.List[str], insertions: typing.List[str], deletions: typing.List[str]) -> typing.List[SeqRecord]:
    if pdb_id is None:
        pdb_id = ""
    if description is None:
        description = ""

    seqs = extract_sequence_from_pdb_file(path)
    structure = parse_to_structure(path)

    if mutations:
        seqs = mutate_residues(structure, seqs, parse_mutations(mutations))

    if insertions:
        seqs = add_indels(structure, seqs, list(map(lambda x: dataclasses.replace(x, type=0), parse_indels(insertions))))

    if deletions:
        seqs = add_indels(structure, seqs, list(map(lambda x: dataclasses.replace(x, type=1), parse_indels(deletions))))

    records = [SeqRecord(seq,
                         id=f'{pdb_id}{idx}' if len(seqs) > 1 else pdb_id,
                         description=f'{description}') for idx, seq in enumerate(seqs)]
    return records


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

    to_fasta_impl(structure, pdb_id, description, mutations, None, None)


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
