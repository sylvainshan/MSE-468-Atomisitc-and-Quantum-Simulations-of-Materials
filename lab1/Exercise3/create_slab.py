"""
Script to generate a slab structure using pymatgen.

Usage:
    python3 script.py <element> <lattice_constant> <num_layers> <num_vacuum_layers>

Arguments:
    element            : Chemical element symbol (e.g., 'Ag')
    lattice_constant   : Lattice constant of the material in angstroms
    num_layers         : Number of atomic layers in the slab
    num_vacuum_layers  : Number of vacuum layers added to the slab

Outputs:
    - LAMMPS data file (e.g., Ag_slab_n_layer4_n_vac3.lmp)
"""

from importlib.util import find_spec
import sys

if find_spec('pymatgen') is None:
    print('pymatgen is not installed.\n')
    print('Please install it by running: pip install --user pymatgen')
    exit()

from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.lammps.data import LammpsData


def create_slab(element: str, alat: float, n_layers: int, n_vac: int):
    """
    Creates an FCC metal slab structure using pymatgen and writes a LAMMPS file.

    Parameters:
        element (str): Chemical element symbol (e.g., 'Ag')
        alat (float): Lattice constant in angstroms
        n_layers (int): Number of atomic layers in the slab
        n_vac (int): Vacuum thickness in angstroms
    """
    # Create conventional FCC bulk structure
    lattice = Lattice.cubic(alat)
    species = [element] * 4
    frac_coords = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
    bulk_structure = Structure(lattice, species, frac_coords)

    # Generate the slab
    slab_gen = SlabGenerator(
        bulk_structure, 
        miller_index=[1, 0, 0], 
        min_slab_size=n_layers,  # Slab thickness
        min_vacuum_size=n_vac,  # Vacuum thickness
        in_unit_planes=True, # Number of layers/vacuum in units of atomic planes instead of Angstroms
        center_slab=False,  # Centers the slab
        primitive=False
    )
    slab = slab_gen.get_slab()

    print("Lattice vectors:\n", slab.lattice)
    print("Atomic positions (fractional):\n", slab.frac_coords)

    # Save to LAMMPS format
    lammps_data = LammpsData.from_structure(slab, atom_style="atomic")
    lammps_data.write_file(f"{element}_slab_nlayer{n_layers}_nvac{n_vac}.lmp")

    print(f"\n{element}_slab_nlayer{n_layers}_nvac{n_vac} has been created.")


if __name__ == '__main__':

    if len(sys.argv) != 5:
        print(
            "Invalid number of arguments.\nCorrect usage: python3 script.py <element> <lattice_constant> <num_layers> <num_vacuum_layers>"
            )
        sys.exit(1)

    element = sys.argv[1]
    alat = float(sys.argv[2])
    n_layers = int(sys.argv[3])
    n_vac = int(sys.argv[4])

    create_slab(element, alat, n_layers, n_vac)