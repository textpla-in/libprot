from ruamel.yaml import YAML

from libprot.pdb import AminoAcid


_YAML = YAML()
_YAML.representer.ignore_aliases = lambda x: True
_YAML.register_class(AminoAcid)
