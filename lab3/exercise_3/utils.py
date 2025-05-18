'''
Utility functions for processing and plotting Density of States (DOS) data.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def load_and_process_nm_dos(filepath, fermi_energy, energy_min, energy_max, skiprows=1):
    """
    Loads and processes non-magnetic (NM) DOS data.

    Args:
        filepath (str): Path to the DOS data file.
        fermi_energy (float): Fermi energy in eV.
        energy_min (float): Minimum energy for cropping (relative to Fermi energy).
        energy_max (float): Maximum energy for cropping (relative to Fermi energy).
        skiprows (int, optional): Number of rows to skip at the beginning of the file. Defaults to 1.

    Returns:
        tuple: (nm_energy_cropped, nm_dos_cropped)
               - nm_energy_cropped (np.ndarray): Cropped energy values (eV), shifted by Fermi energy.
               - nm_dos_cropped (np.ndarray): Cropped DOS values (states/eV).
    """
    data = np.loadtxt(filepath, skiprows=skiprows)
    energy = data[:, 0] - fermi_energy
    dos = data[:, 1]

    mask = (energy >= energy_min) & (energy <= energy_max)
    return energy[mask], dos[mask]


def load_and_process_fm_dos(filepath, fermi_energy, energy_min, energy_max, skiprows=1):
    """
    Loads and processes ferromagnetic (FM) DOS data, including integrated DOS.

    Args:
        filepath (str): Path to the DOS data file.
        fermi_energy (float): Fermi energy in eV.
        energy_min (float): Minimum energy for cropping (relative to Fermi energy).
        energy_max (float): Maximum energy for cropping (relative to Fermi energy).
        skiprows (int, optional): Number of rows to skip at the beginning of the file. Defaults to 1.

    Returns:
        tuple: (fm_energy_cropped, fm_dos_up_cropped, fm_dos_down_cropped, int_dos_up, int_dos_down)
               - fm_energy_cropped (np.ndarray): Cropped energy values (eV), shifted by Fermi energy.
               - fm_dos_up_cropped (np.ndarray): Cropped spin-up DOS values (states/eV).
               - fm_dos_down_cropped (np.ndarray): Cropped spin-down DOS values (states/eV), made negative for plotting.
               - int_dos_up (np.ndarray): Integrated spin-up DOS, renormalized to start at 0.
               - int_dos_down (np.ndarray): Integrated spin-down DOS, renormalized to start at 0.
    """
    data = np.loadtxt(filepath, skiprows=skiprows)
    energy = data[:, 0] - fermi_energy
    dos_up = data[:, 1]
    dos_down = -data[:, 2]  # Make spin down negative for plotting convention

    mask = (energy >= energy_min) & (energy <= energy_max)
    energy_cropped = energy[mask]
    dos_up_cropped = dos_up[mask]
    dos_down_cropped = dos_down[mask]

    # Calculate integrated DOS from cropped data
    delta_e = energy_cropped[1] - energy_cropped[0] if len(energy_cropped) > 1 else 0
    int_dos_up = np.cumsum(dos_up_cropped) * delta_e
    int_dos_down = np.cumsum(dos_down_cropped) * delta_e # dos_down_cropped is already negative

    # Renormalize to start at 0
    int_dos_up = int_dos_up - int_dos_up[0] if len(int_dos_up) > 0 else np.array([])
    int_dos_down = int_dos_down - int_dos_down[0] if len(int_dos_down) > 0 else np.array([])

    return energy_cropped, dos_up_cropped, dos_down_cropped, int_dos_up, int_dos_down


def plot_combined_dos(nm_energy, nm_dos, fm_energy, fm_dos_up, fm_dos_down, int_dos_up, int_dos_down, energy_min_plot, energy_max_plot, nm_ylim_max=10, fm_ylim_abs_max=5.5):
    """
    Plots combined non-magnetic and ferromagnetic DOS, including integrated DOS for FM.

    Args:
        nm_energy (np.ndarray): Energy for NM DOS.
        nm_dos (np.ndarray): NM DOS values.
        fm_energy (np.ndarray): Energy for FM DOS.
        fm_dos_up (np.ndarray): Spin-up FM DOS values.
        fm_dos_down (np.ndarray): Spin-down FM DOS values (should be negative).
        int_dos_up (np.ndarray): Integrated spin-up FM DOS.
        int_dos_down (np.ndarray): Integrated spin-down FM DOS.
        energy_min_plot (float): Minimum energy for x-axis plot limit.
        energy_max_plot (float): Maximum energy for x-axis plot limit.
        nm_ylim_max (float, optional): Max y-limit for NM DOS plot. Defaults to 10.
        fm_ylim_abs_max (float, optional): Absolute max y-limit for FM DOS plot. Defaults to 5.5.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True)

    # Plot 1: Non-magnetic DOS
    ax1.plot(nm_energy, nm_dos, 'b-', label='NM')
    ax1.axvline(x=0, color='k', linestyle='--', label=r'$E_F$')
    ax1.set_ylabel('DOS [states/eV]')
    ax1.set_title('Non-Magnetic HCP Cobalt')
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax1.set_ylim(0, nm_ylim_max)
    ax1.set_xlabel('Energy [eV]')

    # Plot 2: Ferromagnetic DOS
    ax2.fill_between(fm_energy, fm_dos_up, color='coral', alpha=0.2)
    ax2.plot(fm_energy, fm_dos_up, color='red', linestyle='--', label='Spin ↑', alpha=0.8, linewidth=1)
    ax2.fill_between(fm_energy, fm_dos_down, color='skyblue', alpha=0.2) # fm_dos_down is already negative
    ax2.plot(fm_energy, fm_dos_down, 'b--', label='Spin ↓', alpha=0.8, linewidth=1)
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    ax2.set_ylim(-fm_ylim_abs_max, fm_ylim_abs_max)
    yticks = ax2.get_yticks()
    yticklabels = [str(abs(int(y))) for y in yticks]
    ax2.set_yticklabels(yticklabels)
    ax2.set_ylabel('DOS [states/eV]')
    ax2.set_title('Ferromagnetic HCP Cobalt')
    ax2.grid(alpha=0.3)
    ax2.set_xlim(energy_min_plot, energy_max_plot)

    # Create a secondary y-axis for integrated DOS
    ax2_twin = ax2.twinx()
    ax2_twin.plot(fm_energy, int_dos_up, 'r-', label='Int. DOS ↑')
    ax2_twin.plot(fm_energy, int_dos_down, 'b-', label='Int. DOS ↓') # int_dos_down is already negative or starts at 0
    ax2_twin.set_ylabel('Integrated DOS [states]')

    # Adjust y-limits for integrated DOS, ensuring symmetry around 0 if possible
    max_abs_int_dos = 0
    if len(int_dos_up) > 0:
        max_abs_int_dos = max(max_abs_int_dos, np.abs(int_dos_up).max())
    if len(int_dos_down) > 0:
        max_abs_int_dos = max(max_abs_int_dos, np.abs(int_dos_down).max())
    
    if max_abs_int_dos == 0: # Handle case with no integrated DOS data
        max_abs_int_dos = 1 # Default to avoid division by zero or empty plot

    ax2_twin.set_ylim(-max_abs_int_dos, max_abs_int_dos)
    tick_positions = np.linspace(-max_abs_int_dos, max_abs_int_dos, 9)
    ax2_twin.set_yticks(tick_positions)
    y_tick_labels_twin = []
    for y_val in tick_positions:
        if abs(y_val) < 1e-3: # Check for values very close to zero
            y_tick_labels_twin.append("0")
        else:
            y_tick_labels_twin.append(str(abs(int(round(y_val)))))
    ax2_twin.set_yticklabels(y_tick_labels_twin)

    ax2.axvline(x=0, color='k', linestyle='--', alpha=0.5)

    # Create combined legend
    lines_1, labels_1 = ax2.get_legend_handles_labels()
    lines_2, labels_2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')

    # Add a marker at the Fermi energy for both integrated DOS if data exists
    if len(fm_energy) > 0:
        fermi_idx = np.abs(fm_energy).argmin()  # Index of the Fermi energy
        if fermi_idx < len(int_dos_up):
            ax2_twin.scatter(0, int_dos_up[fermi_idx], color='r', marker='o', s=100, 
                            label='Int. DOS at EF ↑', zorder=5, edgecolor='black')
        if fermi_idx < len(int_dos_down):
            ax2_twin.scatter(0, int_dos_down[fermi_idx], color='b', marker='o', s=100, 
                            label='Int. DOS at EF ↓', zorder=5, edgecolor='black')

    ax2.set_xlabel('Energy [eV]')
    plt.tight_layout()
    plt.show()


def set_latex_fonts():
    matplotlib.rcdefaults()
    matplotlib.rcParams['text.usetex'] = True            
    matplotlib.rcParams.update({
        'font.family': 'serif',                         
        'font.serif': 'Palatino',                       
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{mathpazo}',  
    })

    # Bolded frame
    matplotlib.rcParams['axes.labelweight'] = 'bold'  # Axis label weight
    matplotlib.rcParams['axes.titleweight'] = 'bold'  # Axis title weight
    matplotlib.rcParams['axes.linewidth'] = 1.5      # Axis line width
    
    # Fontsize
    matplotlib.rcParams['axes.labelsize'] = 24       # Axis label font size
    matplotlib.rcParams['axes.titlesize'] = 24       # Axis title font size
    matplotlib.rcParams['xtick.labelsize'] = 20      # X-axis tick label size
    matplotlib.rcParams['ytick.labelsize'] = 20      # Y-axis tick label size
    matplotlib.rcParams['legend.fontsize'] = 22      # Legend font size

    # Customize ticks
    matplotlib.rcParams['xtick.major.width'] = 1.2   # X-axis major tick width
    matplotlib.rcParams['ytick.major.width'] = 1.2   # Y-axis major tick width
    matplotlib.rcParams['xtick.minor.size'] = 4      # X-axis minor tick size
    matplotlib.rcParams['ytick.minor.size'] = 4      # Y-axis minor tick size
    matplotlib.rcParams['xtick.major.size'] = 8      # X-axis major tick size
    matplotlib.rcParams['ytick.major.size'] = 8      # Y-axis major tick size
    matplotlib.rcParams['xtick.minor.width'] = 1.2   # X-axis minor tick width
    matplotlib.rcParams['ytick.minor.width'] = 1.2   # Y-axis minor tick width

    # Customize legend
    matplotlib.rcParams['legend.frameon'] = True     # Enable/Disable the frame around the legend
    matplotlib.rcParams['legend.framealpha'] = 1     # Legend frame transparency
    matplotlib.rcParams['legend.loc'] = 'best'       # Legend location
    matplotlib.rcParams['legend.handlelength'] = 1   # Legend handle length
    matplotlib.rcParams['legend.handleheight'] = 1   # Legend handle height
    matplotlib.rcParams['legend.borderpad'] = 0.5    # Legend border padding

    # Customize grid
    matplotlib.rcParams['grid.color'] = 'gray'       # Grid color
    matplotlib.rcParams['grid.linestyle'] = '-'      # Grid line style
    matplotlib.rcParams['grid.linewidth'] = 1      # Grid line width
    matplotlib.rcParams['grid.alpha'] = 0.3          # Grid line transparency