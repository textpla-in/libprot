from __future__ import annotations

import typing
from dataclasses import dataclass, field, asdict
from io import StringIO

from ruamel.yaml import YAML

from ..pdb import ResidueModifier


@dataclass
class Stability:

    osprey_version: float = 3.1
    design_name: str = ""
    epsilon: float = 0.99
    residue_configurations: typing.List[ResidueModifier] = field(default_factory=list)
    molecule: str = ""

    def to_yaml(self) -> str:
        dict_repr = asdict(self)
        yaml = YAML()
        str_out = StringIO()
        yaml.dump(dict_repr, str_out)
        return str_out.getvalue()

    @classmethod
    def from_yaml(cls, stream: typing.TextIO) -> Stability:
        d = YAML().load(stream)
        stability = cls(**d)
        return stability
