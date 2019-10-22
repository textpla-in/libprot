from dataclasses import dataclass
from functools import total_ordering


@dataclass
class Mutation:
    chain: str = ''
    index: int = 0
    from_aa: str = ''
    to_aa: str = ''


@dataclass
@total_ordering
class Indel:
    chain: str = ''
    index: int = 0
    aa: str = ''
    type: int = 0 # insertion = 0, deletion=1

    def is_insertion(self):
        return self.type == 0

    def __lt__(self, other):
        return self.type < other.type or \
               self.chain < other.chain or \
               self.index < other.index or \
               self.aa < other.aa
