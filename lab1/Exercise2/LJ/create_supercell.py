"""
Script to generate a bulk supercell using the Atomic Simulation Environment (ASE).

Usage:
    python3 create_supercell.py <element> <lattice_constant> <nx> <ny> <nz>

Arguments:
    element          : Chemical element symbol (e.g., 'Ag' for Silver)
    lattice_constant : Lattice constant of the material in angstroms
    nx               : Repetition factor along x-axis
    ny               : Repetition factor along y-axis
    nz               : Repetition factor along z-axis

Output:
    - LAMMPS data file (e.g., Ag_supercell_2x2x3.lmp)
"""

from importlib.util import find_spec
import sys

if find_spec('ase') is None:
    print('Atomic Simulation Environment (ASE) is not installed.\n')
    print('Please install it by running: pip install --user ase')
    exit()

from ase.build import bulk
from ase.io.lammpsdata import write_lammps_data

def create_supercell(element: str, alat: float, nx: int, ny: int, nz: int):
    """
    Creates a bulk supercell structure and writes a LAMMPS data file.

    Parameters:
        element (str): Chemical element symbol (e.g., 'Ag')
        alat (float): Lattice constant in angstroms
        nx (int): Repetition factor along x-axis
        ny (int): Repetition factor along y-axis
        nz (int): Repetition factor along z-axis
    """
    bulk_structure = bulk(element, 'fcc', a=alat, cubic=True)
    
    supercell = bulk_structure.repeat((nx, ny, nz))
    
    output_filename = f'{element}_supercell_{nx}x{ny}x{nz}.lmp'
    write_lammps_data(output_filename, supercell, masses=True)
    
    print(f'\n{output_filename} has been created.')

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Invalid number of arguments.\n")
        print("Correct usage: python3 create_supercell.py <element> <lattice_constant> <nx> <ny> <nz>")
        sys.exit(1)
    
    element = sys.argv[1]
    try:
        alat = float(sys.argv[2])
        nx = int(sys.argv[3])
        ny = int(sys.argv[4])
        nz = int(sys.argv[5])
    except ValueError:
        print("Error: 'lattice_constant' must be a float and scaling factors 'nx', 'ny', 'nz' must be integers.")
        sys.exit(1)

    create_supercell(element, alat, nx, ny, nz)