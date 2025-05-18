import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import lmfit
from IPython.display import display, HTML


BOHR = 0.529177210544  # Bohr radius in Angstroms
RY = 13.605693122990   # Rydberg energy in eV


conduction_band_color = '#067eba' 
valence_band_color = '#a146c2'


def plot_single_layer_band_structure(celldm_3=30, fermi_level_offset=-5.0914, num_bands_to_plot=60, show=True): 

    a = 5.97 # Bohr
    c = a * celldm_3
    print("="*50)
    print(f"Plotting band structure for celldm(3) = c/a = {celldm_3}")
    print(f"This corresponds to a vacuum thickness of {c:.4f} Bohr")
    print("="*50)

    # Load the data
    filepath = "single_layer/celldm(3) variation/"
    data = np.loadtxt(filepath + "MoS2_" + str(celldm_3) + ".dat.gnu")
    k = np.unique(data[:, 0])
    bands = np.reshape(data[:, 1], (-1, len(k)))

    # Define the k-path
    Gamma = np.array([0, 0, 0])
    K = np.array([1/3, 1/3, 0])
    M = np.array([0.5, 0, 0])
    K_to_Gamma = np.linalg.norm(K - Gamma)
    Gamma_to_M = np.linalg.norm(Gamma - M)
    M_to_K = np.linalg.norm(K - M)
    
    k_points = [
        0, 
        K_to_Gamma,
        K_to_Gamma + Gamma_to_M,
        K_to_Gamma + Gamma_to_M + M_to_K
    ]

    labels = [
        r"K",
        r"$\Gamma$",
        r"M",
        r"K"
    ]

    # Find the valence band maxima and conduction band minima
    valence_band_index = 12
    conduction_band_index = 13
    vbm = np.max(bands[valence_band_index, :])
    cbm = np.min(bands[conduction_band_index, :])
    band_gap = cbm - vbm
    print(f"Overall VBM: {vbm - fermi_level_offset:.4f} eV")
    print(f"Overall CBM: {cbm - fermi_level_offset:.4f} eV")
    print(f"Overall Band gap: {band_gap:.4f} eV")

    vbm_idx = np.argmax(bands[valence_band_index, :])
    cbm_idx = np.argmin(bands[conduction_band_index, :])
    vbm_k = k[vbm_idx]
    cbm_k = k[cbm_idx]
    
    plt.figure(figsize=(8, 6))
    for band in range(len(bands)):
        if band < num_bands_to_plot:
            plt.plot(k, bands[band, :] - fermi_level_offset,
                    linewidth=1.5,
                    alpha=0.8,
                    color=valence_band_color if band <= valence_band_index else conduction_band_color)

    plt.scatter(vbm_k, vbm - fermi_level_offset, color=valence_band_color, s=100, zorder=5, marker='o', label='VBM')
    plt.scatter(cbm_k, cbm - fermi_level_offset, color=conduction_band_color, s=100, zorder=5, marker='o', label='CBM')
                    
    # plt.axhline(y=0, color='k', linestyle='--', linewidth=1.5, alpha=0.5, label='Fermi level')
    plt.xlabel(r"$\textbf{Wavevector}$")
    plt.xticks(k_points, labels)
    plt.xlim(0, k_points[-1])
    plt.ylabel(r"$\textbf{Energy [eV]}$")
    plt.grid(axis='x', linestyle='-', linewidth=1, alpha=0.7, color='k')
    plt.grid(axis='y', linestyle='', linewidth=1, alpha=0.5)
    plt.ylim(-6, 6)
    plt.tight_layout()

    if show:
        plt.show()

    print("="*50)

    return band_gap


def parabola(x, amplitude, center, offset):
    """Defines a parabola: f(x) = amplitude * (x - center)**2 + offset"""
    return amplitude * (x - center)**2 + offset

def fit_parabola(lattice_parameter, energy): 
    parabolic_model = lmfit.Model(parabola)

    min_idx = np.argmin(energy)
    initial_center = lattice_parameter[min_idx]
    initial_offset = energy[min_idx]
    initial_amplitude = (energy[min_idx] - initial_offset) / ((lattice_parameter[0] - initial_center)**2 + 1e-9)
    params = parabolic_model.make_params(amplitude=initial_amplitude if initial_amplitude > 0 else 1,
                                        center=initial_center,
                                        offset=initial_offset)
    result = parabolic_model.fit(energy, params, x=lattice_parameter)
    print(result.fit_report())

    a_fit = result.params['center'].value
    err_a = np.sqrt(result.covar[1, 1])
    E0_fit = result.params['offset'].value
    err_E0 = np.sqrt(result.covar[2, 2])
    amplitude_fit = result.params['amplitude'].value
    err_amplitude = np.sqrt(result.covar[0, 0])

    x_fit = np.linspace(lattice_parameter.min(), lattice_parameter.max(), 1000)
    y_fit = result.eval(x=x_fit)

    display(HTML("<h3>Results</h3>"))
    display(HTML(f"<b>Optimized lattice parameter: ({a_fit:.4f} ± {err_a:.4f}) Å</b>"))
    display(HTML(f"<b>Optimized lattice parameter: ({a_fit / BOHR:.4f} ± {err_a / BOHR:.4f}) Bohr</b>"))
    display(HTML(f"<b>Minimum energy (E_0): ({E0_fit:.4f} ± {err_E0:.4f}) eV</b>"))
    display(HTML(f"<b>Fitted amplitude (A): ({amplitude_fit:.4f} ± {err_amplitude:.4f}) eV/Å²</b>"))

    return x_fit, y_fit



# ======= LATEX =======
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