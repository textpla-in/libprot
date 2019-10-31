import io

from dataclasses import dataclass, field
from enum import Enum

import typing
from prody.proteins import parsePDBStream
from typing import List

from libprot.types import ResidueRenumberingDict


class AminoAcid(Enum):
    ALA = 'Ala'
    ARG = 'Arg'
    ASN = 'Asn'
    ASP = 'Asp'
    CYS = 'Cys'
    GLN = 'Gln'
    GLU = 'Glu'
    GLY = 'Gly'
    HIS = 'His'
    HIE = 'Hie'
    HIP = 'Hip'
    HID = 'Hid'
    ILE = 'Ile'
    LEU = 'Leu'
    LYS = 'Lys'
    MET = 'Met'
    PHE = 'Phe'
    PRO = 'Pro'
    SER = 'Ser'
    THR = 'Thr'
    TRP = 'Trp'
    TYR = 'Tyr'
    VAL = 'Val'
    ASX = 'Asx'
    GLX = 'Glx'
    XAA = 'Xaa'
    TERM = 'TERM'
    OTHER = 'OTHER'

    @classmethod
    def to_yaml(cls, representer, node):
        return representer.represent_str(node.value)

    @classmethod
    def from_yaml(cls, constructor, node):
        return cls[node.value.upper()]


@dataclass(frozen=True)
class Residue:
    chain: str
    res_num: int
    aa_type: AminoAcid
    changed: bool = False


@dataclass(frozen=True)
class Flexibility:
    is_flexible: bool = False
    include_structure_rotamer: bool = False

    def is_default(self):
        return not self.is_flexible and not self.include_structure_rotamer


@dataclass
class ResidueModifier:
    identity: Residue
    flexibility: Flexibility
    mutability: typing.List[AminoAcid] = field(default_factory=list)


def get_amino_acids(stream: typing.TextIO) -> List[Residue]:
    structure = parsePDBStream(stream).getHierView()

    rv = []
    for x in structure.iterResidues():
        try:
            aa_type = AminoAcid[x.getResname().upper()]
        except KeyError:
            aa_type = AminoAcid.OTHER

        rv.append(Residue(chain=str(x.getChid()), res_num=int(x.getResnum()), aa_type=aa_type))

    return rv


def pad_line(line):
    """Helper function to pad line to 80 characters in case it is shorter"""
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]  # 80 + newline character


def renumber_pdb(pdb: typing.TextIO, d: ResidueRenumberingDict) -> typing.TextIO:
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')

    lines = []
    for line in [pad_line(line) for line in pdb]:
        if line.startswith(records):
            replacement = d[int(line[22:26])]
            lines.append(line[:22] + str(replacement).rjust(4) + line[26:])
        else:
            lines.append(line)
    return io.StringIO(''.join(lines))


