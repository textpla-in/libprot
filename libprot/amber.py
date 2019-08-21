import io
import os
import subprocess
from enum import Enum
from pathlib import Path

_amber_home = ''


def ensure_ambertools():
    # AMBERHOME should be defined in the environment
    global _amber_home
    env = os.environ
    amberhome = "AMBERHOME"
    if not env.get(amberhome):
        raise EnvironmentError(f"Expecting {amberhome} to be defined in the user's environment")
    _amber_home = env[amberhome]


class ForceFieldType(Enum):
    PROTEIN = 'leaprc.protein.ff14SB'
    RNA = 'leaprc.RNA.OL3'
    DNA = 'leaprc.DNA.OL15'
    WATER = 'leaprc.water.tip3p'  # not needed if using implicit solvent model
    GENERAL_FF = 'leaprc.gaff2'


def ensure_bin(bin: Path):
    if not bin.exists():
        raise FileNotFoundError(f"{bin} does not exist")
    if not os.access(str(bin), os.X_OK):
        raise PermissionError(f"{bin} is not executable")


def run_subprocess(bin: Path, stdin, args):
    ensure_bin(bin)
    completed_proc = subprocess.run([bin] + args, stdin=stdin, capture_output=True, text=True, check=True)
    return completed_proc.returncode, io.StringIO(completed_proc.stdout), io.StringIO(completed_proc.stderr)


def prep_pdb_for_amber(instream: io.TextIOBase):
    """Runs pdb4amber on the structure
    :param instream: an instream for the PDB.
    :return a readable stream for the prepped PDB
    """
    path_to_bin = Path(_amber_home).joinpath("bin", "pdb4amber")
    return_code, outstream, errstream = run_subprocess(path_to_bin, instream, ["--dry"])
    return outstream


ensure_ambertools()
