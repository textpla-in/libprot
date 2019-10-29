from itertools import chain
from pathlib import Path

from libprot import cli
from libprot.cli import mutate_residues, extract_seqs_from_pdb_file, parse_to_structure, \
    add_indels, mutate_and_indels
from libprot.types import Indel, Mutation


def test_parses_insertions():
    insertions = ['A:S172', 'A:M255']
    parsed = cli.parse_indels(insertions)
    assert parsed == [Indel('A', 172, 'S'), Indel('A', 255, 'M')]


def test_parses_mutations():
    mutations = ['A:T79M', 'B:S101A']
    parsed = cli.parse_mutations(mutations)
    assert parsed == [Mutation('A', 79, 'T', 'M'), Mutation('B', 101, 'S', 'A')]


def test_seq_from_structure():
    sequences = extract_seqs_from_pdb_file(Path('resources/6ODG.pdb'))
    assert len(sequences) == 2
    assert (len(list(chain.from_iterable(sequences)))) == 12


def test_mutate():
    mut1 = Mutation('A', 1, 'S', 'M')
    mut2 = Mutation('B', 1, 'S', 'M')
    structure_path = Path('resources/6ODG.pdb')
    seqs = extract_seqs_from_pdb_file(structure_path)
    new_seqs = mutate_residues(parse_to_structure(structure_path), seqs, [mut1, mut2])
    assert new_seqs[0][0] == 'M'
    assert new_seqs[1][0] == 'M'


def test_add_insertion():
    ins1 = Indel('A', 1, 'M', type=0)
    ins2 = Indel('B', 1, 'M', type=0)
    ins3 = Indel('A', 6, 'M', type=0)
    structure_path = Path('resources/6ODG.pdb')
    seqs = extract_seqs_from_pdb_file(structure_path)
    new_seqs = add_indels(parse_to_structure(structure_path), seqs, [ins1, ins2, ins3])
    assert len(new_seqs[0]) == 8
    assert new_seqs[0][0] == 'M'
    assert new_seqs[0][5] == 'M'
    assert len(new_seqs[1]) == 7
    assert new_seqs[1][0] == 'M'


def test_mutate_and_insert():
    ins1 = Indel('A', 1, 'G', type=0)
    mut1 = Mutation('A', 1, 'S', 'M')
    structure_path = Path('resources/6ODG.pdb')
    seqs = extract_seqs_from_pdb_file(structure_path)
    new_seqs = mutate_and_indels(parse_to_structure(structure_path), seqs, [mut1], [ins1])
    assert len(new_seqs) == 2
    assert len(new_seqs[0]) == 7
    for i, aa in enumerate('GMVQIVY'):
        assert new_seqs[0][i] == aa
