import operator
import typing
from dataclasses import dataclass
from pathlib import Path
import functools

from Bio.PDB import PDBParser, Structure
from prody import parsePDB, AtomGroup


@dataclass(frozen=True)
class ResidueIdentifier:
    res_name: str
    res_seq: int


def find_mutated_residues(template: Structure, new_struct: Structure) -> typing.Set[ResidueIdentifier]:
    generator = zip(template.get_residues(), new_struct.get_residues())
    return {ResidueIdentifier(struct_res.get_resname(), int(struct_res.get_id()[1]))
            for temp_res, struct_res in generator
            if temp_res.get_resname() != struct_res.get_resname()}


def get_residues_in_shell(pdb_path: Path, residue: ResidueIdentifier, shell_angs: float) -> typing.Set[ResidueIdentifier]:
    atom_group = parsePDB(pdb_path)
    return _get_residues_in_shell(atom_group, residue, shell_angs)


def _get_residues_in_shell(atom_group: AtomGroup, residue: ResidueIdentifier, shell_angs: float) -> typing.Set[ResidueIdentifier]:
    residue_atoms = atom_group.select(f'resnum {residue.res_seq}')
    shell_atoms = atom_group.select(f'same residue as exwithin {shell_angs} of ag', ag=residue_atoms)
    if not shell_atoms:
        return set()

    return set(ResidueIdentifier(res_name=atom.getResname(), res_seq=int(atom.getResnum())) for atom in shell_atoms)


def mutants_and_residues_in_shell(template_structure: Path, target_structure: Path, shell: float = 0.0) -> typing.Set[ResidueIdentifier]:

    with template_structure.open() as template, target_structure.open() as target:
        # Find mutations
        bio_template = PDBParser(QUIET=True).get_structure(template_structure.name, template)
        bio_target = PDBParser(QUIET=True).get_structure(target_structure.name, target)
        diffs = find_mutated_residues(bio_template, bio_target)

        if shell == 0.0:
            return diffs

        # Find all residues within any of the mutations
        prody_target = parsePDB(target_structure)
        in_shell_residues = [_get_residues_in_shell(prody_target, flexible_residues, shell)
                             for flexible_residues in diffs]
        all_in_shells = functools.reduce(operator.ior, in_shell_residues)

        return all_in_shells | diffs
