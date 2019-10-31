import io
import os
import subprocess
import sys
import tempfile

import typing
from enum import Enum
from pathlib import Path

from libprot.pdb import renumber_pdb
from libprot.types import ResidueRenumberingDict


## This generates the topology and coordinates using the protein force field
# source leaprc.protein.ff14SB
# pdb = loadPdb /home/nsg/tmp/KE07_apo_pdb4amber.pdb
# saveamberparm pdb prmtop prmcrd
# quit


# See http://ambermd.org/tutorials/basic/tutorial1/section3.htm for invocation and options of sander
## This runs sander
# sander -O -i mdin -o mdout -c rst7 -p prmtop -r ncrst
## mdin:
# Test run 1
#  &cntrl
#     IMIN = 1,
#     NCYC = 250,
#     MAXCYC = 500,
#     NTB = 0,
#     IGB = 0,
#     CUT = 12
# /

class ForceFieldType(Enum):
    PROTEIN = 'leaprc.protein.ff14SB'
    RNA = 'leaprc.RNA.OL3'
    DNA = 'leaprc.DNA.OL15'
    WATER = 'leaprc.water.tip3p'  # not needed if using implicit solvent model
    GENERAL_FF = 'leaprc.gaff2'


def make_amber_env():
    """AMBERHOME must be a valid environment variable, otherwise this function will throw.
    :returns an Amber environment, which is copied from amber.sh in the amber distribution"""
    amberhome = "AMBERHOME"
    env = dict(os.environ)
    if not env.get(amberhome):
        raise EnvironmentError(f"Expecting {amberhome} to be defined in the user's environment")

    env['AMBER_PREFIX'] = amber_prefix = env[amberhome]
    env['PATH'] = f"{amber_prefix}:{env.get('PATH')}"
    env['PYTHONPATH'] = f"{amber_prefix}/lib/python2.7/site-packages:{env.get('PYTHONPATH')}"
    env['LD_LIBRARY_PATH'] = f"{amber_prefix}/lib:{env.get('LD_LIBRARY_PATH')}"
    return env


def ensure_bin(bin: Path):
    if not bin.exists():
        raise FileNotFoundError(f"{bin} does not exist")
    if not os.access(str(bin), os.X_OK):
        raise PermissionError(f"{bin} is not executable")


def run_subprocess(bin: Path, stdin, args) -> typing.Tuple[int, io.StringIO, io.StringIO]:
    completed_proc = subprocess.run([str(bin)] + args, stdin=(stdin or sys.stdin), capture_output=True, text=True,
                                    env=make_amber_env())
    return completed_proc.returncode, io.StringIO(completed_proc.stdout), io.StringIO(completed_proc.stderr)


def get_amber_bin(program_name: str) -> Path:
    """Finds a program in the bin directory """
    env = make_amber_env()
    executable = Path(env['AMBERHOME']).joinpath("bin", program_name)
    ensure_bin(executable)
    return executable


def parse_renumbering_file(name: str) -> ResidueRenumberingDict:
    with open(name) as f:

        renumberings = {}
        for line in f:
            from_aa, from_num, to_aa, to_num = line.split()
            renumberings[int(to_num)] = int(from_num)

        return renumberings


def prep_pdb_for_amber(instream: typing.TextIO) -> typing.Tuple[int, typing.TextIO, typing.TextIO, ResidueRenumberingDict]:
    """Runs pdb4amber on the structure
    :param instream: an instream for the PDB.
    :return tuple(return code, stdout, stderr)
    """
    executable = get_amber_bin("pdb4amber")
    ret_code, stdout, stderr = run_subprocess(executable, instream, ["--dry"])
    renumbers = parse_renumbering_file('stdout_renum.txt')
    return ret_code, stdout, stderr, renumbers


def reduce_pdb(instream: typing.TextIO, remove_hydrogens: bool = False) -> typing.Tuple[int, typing.TextIO, typing.TextIO]:
    arg = '-trim' if remove_hydrogens else '-build'
    executable = get_amber_bin("reduce")
    return run_subprocess(executable, instream, [arg, '-'])


def make_parameter_files(in_stream: typing.TextIO, ff_type: ForceFieldType) -> typing.Tuple[
    int, typing.TextIO, typing.TextIO]:
    # write input structure to temporary file
    tmp_pdb_fd, tmp_pdb_path = tempfile.mkstemp(text=True)
    tmp_prmtop_fd, tmp_prmtop_path = tempfile.mkstemp(text=True)
    tmp_prmcrd_fd, tmp_prmcrd_path = tempfile.mkstemp(text=True)
    tmp_tleap_fd, tmp_tleap_path = tempfile.mkstemp(text=True)

    tleap_input = f"""source {ff_type.value}
pdb = loadPdb {tmp_pdb_path}
saveamberparm pdb {tmp_prmtop_path} {tmp_prmcrd_path}
quit"""

    with open(tmp_pdb_fd, mode='r+') as tmp_pdb, \
        open(tmp_prmtop_fd, mode='r+') as tmp_prmtop, \
        open(tmp_prmcrd_fd, mode='r+') as tmp_prmcrd, \
        open(tmp_tleap_fd, mode='r+') as tmp_tleap:
        tmp_tleap.write(tleap_input)
        tmp_tleap.seek(0)
        tmp_pdb.write(in_stream.read())

        executable = get_amber_bin("tleap")
        return_code, stdout, stderr = run_subprocess(executable, tmp_tleap, ["-s", "-f", "-"])
        tmp_prmtop.seek(0)
        tmp_prmcrd.seek(0)
        prmtop = io.StringIO(tmp_prmtop.read())
        prmcrd = io.StringIO(tmp_prmcrd.read())

    for path in [tmp_pdb_path, tmp_prmtop_path, tmp_prmcrd_path, tmp_tleap_path]:
        os.unlink(path)

    return return_code, prmtop, prmcrd


def run_sander(topology: typing.TextIO, coordinates: typing.TextIO, max_cycles: int) -> typing.Tuple[
    int, typing.BinaryIO, typing.TextIO]:
    """Runs sander on the PDB to minimize it"""
    tmp_topology_fd, tmp_topology_path = tempfile.mkstemp(text=True)
    tmp_coordinates_fd, tmp_coordinates_path = tempfile.mkstemp(text=True)
    tmp_mdin_fd, tmp_mdin_path = tempfile.mkstemp(text=True)
    tmp_restart_fd, tmp_restart_path = tempfile.mkstemp(text=True)
    tmp_mdout_fd, tmp_mdout_path = tempfile.mkstemp(text=True)

    input_file = f"""
Minimization with implicit solvent
 &cntrl
    IMIN = 1,
    MAXCYC = {max_cycles},
    IGB = 1,
    CUT = 12
/
"""

    with open(tmp_topology_fd, 'w') as tmp_topology, open(tmp_coordinates_fd, 'w') as tmp_coordinates, open(tmp_mdin_fd,
                                                                                                            'w') as tmp_mdin:
        tmp_topology.write(topology.read())
        tmp_coordinates.write(coordinates.read())
        tmp_mdin.write(input_file)

    executable = get_amber_bin("sander")
    args = ['-O', '-i', tmp_mdin_path, '-o', tmp_mdout_path, '-c', tmp_coordinates_path, '-p', tmp_topology_path, '-r',
            tmp_restart_path]

    try:
        return_code, stdout, stderr = run_subprocess(executable, None, args)
        with open(tmp_restart_fd, 'rb') as final_coords:
            return return_code, io.BytesIO(final_coords.read()), stderr
    finally:
        for path in [tmp_topology_path, tmp_coordinates_path, tmp_mdin_path, tmp_restart_path, tmp_mdout_path]:
            os.unlink(path)


def convert_coords_to_pdb(topology: typing.TextIO, coordinates: typing.BinaryIO) -> typing.Tuple[int, typing.TextIO, typing.TextIO]:
    tmp_topology_fd, tmp_topology_path = tempfile.mkstemp(text=True)
    tmp_coordinates_fd, tmp_coordinates_path = tempfile.mkstemp(text=True)
    tmp_pdb_fd, tmp_pdb_path = tempfile.mkstemp(text=True)

    with open(tmp_topology_fd, 'w') as tmp_topology, open(tmp_coordinates_fd, 'wb') as tmp_coordinates:
        tmp_topology.write(topology.read())
        tmp_coordinates.write(coordinates.read())

    executable = get_amber_bin("ambpdb")

    try:
        return run_subprocess(executable, None, ["-p", tmp_topology_path, "-c", tmp_coordinates_path])
    finally:
        for path in [tmp_topology_path, tmp_coordinates_path, tmp_pdb_path]:
            os.unlink(path)


def do_minimize_pdb(pdb: typing.TextIO, ff_type: ForceFieldType, rounds: int) -> typing.TextIO:
    return_code, stdout, stderr, renumber_dict = prep_pdb_for_amber(pdb)
    return_code, prmtop, prmcrd = make_parameter_files(stdout, ff_type)
    topology = prmtop.read()
    return_code, final_coords, stderr = run_sander(io.StringIO(topology), prmcrd, rounds)
    return_code, new_pdb, stderr = convert_coords_to_pdb(io.StringIO(topology), final_coords)
    return renumber_pdb(new_pdb, renumber_dict)
