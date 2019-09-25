from __future__ import annotations

import typing
from dataclasses import dataclass, field, asdict
from io import StringIO

from ruamel.yaml.scalarstring import PreservedScalarString as L

from ..pdb import ResidueModifier

from . import _YAML


@dataclass
class Stability:

    osprey_version: float = 3.1
    design_name: str = ""
    epsilon: float = 0.99
    residue_configurations: typing.List[ResidueModifier] = field(default_factory=list)
    molecule: str = ""

    def to_yaml(self) -> str:
        dict_repr = asdict(self)
        dict_repr["molecule"] = L(self.molecule)
        str_out = StringIO()
        _YAML.dump(dict_repr, str_out)
        return str_out.getvalue()

    @classmethod
    def from_yaml(cls, stream: typing.TextIO) -> Stability:
        d = _YAML.load(stream)
        stability = cls(**d)
        return stability
