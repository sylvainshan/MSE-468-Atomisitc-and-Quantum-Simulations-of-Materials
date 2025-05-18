import re
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


Ry = 13.605693122990 # eV
Bohr = 0.529177210544 # Angstrom
angstrom = 1e-10  # meter 
eV = 1.602176634e-19  # Joule


def set_mpl_params():
    matplotlib.rcdefaults()
    matplotlib.rcParams['text.usetex'] = True            # Use LaTeX for text rendering

    # Update font settings
    matplotlib.rcParams.update({
        'font.family': 'serif',                          # Use serif font family
        'font.serif': 'Palatino',                        # Use Palatino as the standard font
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{mathpazo}',  # Use the amsmath and mathpazo package for LaTeX
    })

    # Customize the figure size
    matplotlib.rcParams['figure.figsize'] = (8, 6)   # Set the default figure size

    # Customize axes
    matplotlib.rcParams['axes.labelsize'] = 22       # Axis label font size
    matplotlib.rcParams['axes.titlesize'] = 22       # Axis title font size
    matplotlib.rcParams['axes.titlepad'] = 15        # Axis title padding
    matplotlib.rcParams['axes.linewidth'] = 1.5      # Axis line width

    # Customize ticks
    matplotlib.rcParams['xtick.labelsize'] = 20      # X-axis tick label size
    matplotlib.rcParams['ytick.labelsize'] = 20      # Y-axis tick label size
    matplotlib.rcParams['xtick.major.width'] = 1.2   # X-axis major tick width
    matplotlib.rcParams['ytick.major.width'] = 1.2   # Y-axis major tick width
    matplotlib.rcParams['xtick.minor.size'] = 4      # X-axis minor tick size
    matplotlib.rcParams['ytick.minor.size'] = 4      # Y-axis minor tick size
    matplotlib.rcParams['xtick.major.size'] = 8      # X-axis major tick size
    matplotlib.rcParams['ytick.major.size'] = 8      # Y-axis major tick size

    # Customize legend
    matplotlib.rcParams['legend.fontsize'] = 20      # Legend font size
    matplotlib.rcParams['legend.frameon'] = True     # Enable/Disable the frame around the legend

    # Customize grid
    matplotlib.rcParams['grid.color'] = 'gray'       # Grid color
    matplotlib.rcParams['grid.linestyle'] = ':'      # Grid line style
    matplotlib.rcParams['grid.linewidth'] = 0.5      # Grid line width

    # Customize lines
    matplotlib.rcParams['lines.linewidth'] = 2.5       # Line width
    matplotlib.rcParams['lines.markersize'] = 10       # Marker size

    # Change figure and axes background colors
    matplotlib.rcParams['figure.facecolor'] = 'white'    # Figure background color
    matplotlib.rcParams['axes.facecolor'] = 'white'      # Axes background color
    

def extract_forces(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    e_cut = None
    atom_data = []
    
    for line in lines:
        if "ecut=" in line:
            match_ecut = re.search(r"ecut=([\d\.]+)", line)
            if match_ecut:
                e_cut = float(match_ecut.group(1))

        if "atom" in line and "force" in line:
            match = re.search(r"atom\s+(\d+)\s+type\s+(\d+)\s+force\s+=\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", line)
            if match:
                atom = int(match.group(1))
                atom_type = int(match.group(2))
                force1 = float(match.group(3))
                force2 = float(match.group(4))
                force3 = float(match.group(5))
                atom_data.append([e_cut, atom, atom_type, force1, force2, force3])
    
    return atom_data


def extract_forces_k(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        
    k = None
    atom_data = []
    for line in lines:
        if "k=" in line:
            match_k = re.search(r"k=([\d\.]+)", line)
            if match_k:
                k = float(match_k.group(1))
        if "atom" in line and "force" in line:
            match = re.search(r"atom\s+(\d+)\s+type\s+(\d+)\s+force\s+=\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", line)
            if match:
                atom = int(match.group(1))
                atom_type = int(match.group(2))
                force1 = float(match.group(3))
                force2 = float(match.group(4))
                force3 = float(match.group(5))
                atom_data.append([k, atom, atom_type, force1, force2, force3])
    return atom_data


def extract_energies(filepath):
    with open(filepath, "r") as f:
        data = f.read()
    pattern = r"ecut=(\d+)\.k=\d+.*?=\s*(-?\d+\.\d+)"
    matches = re.findall(pattern, data)

    df = pd.DataFrame(matches, columns=["E_cut [Ry]", "E [Ry]"])
    df = df.astype({"E_cut [Ry]": int, "E [Ry]": float})
    df_sorted = df.sort_values(by="E_cut [Ry]").reset_index(drop=True)
    df_sorted["E [Ry]"] = df_sorted["E [Ry]"] / 2  # Energy per atom
    return df_sorted


def merge_data(df_1, df_2): 
    df = pd.merge(df_1, df_2, on="E_cut [Ry]", suffixes=("_a9.00", "_a9.05"))

    # Compute energy difference
    df["Delta E [Ry]"] = np.abs(df["E [Ry]_a9.05"] - df["E [Ry]_a9.00"])

    # Conversion to eV
    df["E_cut [eV]"] = df["E_cut [Ry]"] * Ry
    df["Delta E [eV]"] = df["Delta E [Ry]"] * Ry

    # Convergence check
    energy_threshold = 5e-3  # meV/atom
    df["Convergence"] = df["Delta E [eV]"] < energy_threshold

    return df


def extract_energies_a_vary(filepath):
    with open(filepath, "r") as f:
        data = f.read()
    pattern = r"a=(\d+\.\d+).*?total energy\s*=\s*(-?\d+\.\d+)"
    matches = re.findall(pattern, data)
    df = pd.DataFrame(matches, columns=["a [Bohr]", "E [Ry]"])
    df = df.astype({"a [Bohr]": float, "E [Ry]": float})
    df = df.sort_values(by="a [Bohr]").reset_index(drop=True)
    df["E [eV]"] = df["E [Ry]"] * Ry
    df["a [Ang]"] = df["a [Bohr]"] * Bohr
    df["V [Ang^3]"] = df["a [Ang]"] ** 3 / 4 # Volume of the primitive cell 
    return df


def birch_murnaghan(V, E0, B0, B0_prime, V0):
    return E0 + (9/8)*B0*V0*((V0/V)**(2/3) - 1)**2 + (9/16)*B0*V0*(B0_prime - 4)*((V0/V)**(2/3) - 1)**3


def parse_energy_file(filepath):
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            match = re.search(r'scf\.x=([0-9]+\.[0-9]+)(?=\.ecut).*total energy\s*=\s*([-0-9.]+)', line)
            if match:
                x_val = float(match.group(1))
                energy = float(match.group(2))
                data.append({'x': x_val, 'E [Ry]': energy})
        df = pd.DataFrame(data)
        df['iteration'] = df.groupby('x').cumcount() + 1
        df = df.set_index(['x', 'iteration'])
        df = df.sort_index(level=['x', 'iteration'])
    return df

def plot_energy_evolution(df, x_value):
    try:
        df_x = df.loc[x_value]
    except KeyError:
        print(f"No data for x = {x_value}")
        return
    
    plt.figure(figsize=(8, 5))
    plt.plot(df_x.index, df_x['E [Ry]'], marker='o', linestyle='-', color='b')
    plt.title(fr"Energy minimization for $x = {x_value}$")
    plt.xlabel("Iteration number")
    plt.ylabel("Energy [Ry]")
    plt.grid(True)
    plt.show()


def compute_Delta_E(df):
    # Reference energy for the undistorted structure
    E_0 = df.loc[0.00]["E [Ry]"].values[0]

    # For each x value, choose the minimum energy and put them in a new dataframe
    df_E = df.groupby(level=0).min().reset_index()
    df_E["Delta E [Ry]"] = df_E["E [Ry]"] - E_0
    df_E["Delta E [eV]"] = df_E["Delta E [Ry]"] * Ry

    x = df_E["x"].values
    delta_E = df_E["Delta E [eV]"].values

    # The function is symmetric around x=0
    x = np.concatenate((-x[:0:-1], x))
    delta_E = np.concatenate((delta_E[:0:-1], delta_E))
    return x, delta_E


def second_order_polynomial(x, a, b, c):
    return a * x**2 + b * x + c