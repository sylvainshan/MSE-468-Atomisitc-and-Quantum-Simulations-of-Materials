import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import os
import lmfit
from lmfit import report_ci, conf_interval, Model
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from scipy import stats
from scipy.constants import k, atomic_mass
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns
from scipy.constants import Boltzmann, atomic_mass, elementary_charge, gas_constant
import json

# ============================
#     PLOTTING SETTINGS
# =============================

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


# ============================
#     TIMESTEP CONVERGENCE
# =============================

def load_data(filepath): 
    data = pd.read_csv(filepath, sep="\t", comment="#", header=None)
    data.columns = ['t (ps)', 'T (K)', 'E (eV)', 'K (eV)', 'P (Bar)', 'MSD_I', 'MSD_Ag']
    return data

def create_filepaths(timesteps):
    filepaths = {}
    for t in timesteps:
        # Format t with fixed notation 
        t_str = f"{t:.6f}".rstrip('0').rstrip('.') if '.' in f"{t:.6f}" else f"{t:.6f}"
        filepaths[t] = f"P1/timestep_convergence/E_of_t_2_800_{t_str}.dat"
    return filepaths

def compute_total_energy(timesteps, filepaths):
    data = {}
    total_energies = {}
    for t in timesteps:
        try:
            data[t] = load_data(filepaths[t])
            total_energies[t] = data[t]["K (eV)"] + data[t]["E (eV)"]
        except FileNotFoundError:
            print(f"Warning: File {filepaths[t]} not found.")
    return data, total_energies

def plot_std_devs(timesteps, total_energies, convergence_threshold=1e-5, n_atoms=32, verbose=True):
    std_devs = {t: np.std(total_energies[t])/n_atoms for t in timesteps if t in total_energies}
    df = pd.DataFrame.from_dict(std_devs, orient='index', columns=['Standard Deviation'])
    df["Converged"] = df['Standard Deviation'] < convergence_threshold
    
    plt.figure(figsize=(6, 5.3))
    converged = df[df['Converged']]
    not_converged = df[~df['Converged']]

    if verbose: 
        print("Converged timesteps:")
        print(converged)
        print("\nNot converged timesteps:")
        print(not_converged)

    plt.loglog(converged.index, converged['Standard Deviation'], 'o', color="#8b1616", markersize=10, markerfacecolor="#df3113", markeredgewidth=2, zorder=2)
    plt.loglog(not_converged.index, not_converged['Standard Deviation'], 'o', color="#006ead", markersize=10, markerfacecolor='none', markeredgewidth=2, zorder=1)
    plt.loglog(df.index, df['Standard Deviation'], '-', color='black', alpha=0.5, zorder=1)
    plt.xlabel(r"$\boldsymbol{\Delta t}$ $\textbf{[ps]}$")
    plt.ylabel(r"$\boldsymbol{\sigma_E}$ $\textbf{[eV/atom]}$")
    plt.axhline(y=1e-5, color="#000000", linestyle='--', label=rf'Convergence threshold', linewidth=3, zorder=1)
    plt.grid(True)
    plt.legend()

    os.makedirs("Figures", exist_ok=True)
    plt.savefig("Figures/timestep_convergence.pdf", bbox_inches='tight')
    plt.show()


# ============================
#     SUPERCELL CONVERGENCE
# =============================
def create_filepaths_supercell(supercell_sizes, T=283, timestep=0.002):
    filepaths = {}
    for s in supercell_sizes:
        s_str = f"{s:.6f}".rstrip('0').rstrip('.') if '.' in f"{s:.6f}" else f"{s:.6f}"
        filepaths[s] = f"P1/supercell_convergence/T={T}/md_{s_str}_{T}_{timestep}.out"
    return filepaths

def read_supercell_data(filepaths): 
    data = {}
    temperatures = {}
    for s, filepath in filepaths.items():
        try:
            df = pd.read_csv(filepath, sep="\s+", header=None, skiprows=149, nrows=60)
            df.columns = ['Time', 'Temp', 'PotEng', 'KinEng', 'Press', 'MSD_I', 'MSD_Ag']
            data[s] = df

            temperatures[s] = df["Temp"].copy()
        except FileNotFoundError:
            print(f"File not found: {filepath}")
            data[s] = None
        
    return data, temperatures

def get_std_devs(data, supercell_sizes, T=293):
    std_devs = pd.DataFrame.from_dict({s: data[s]["Temp"].std() for s in supercell_sizes if data[s] is not None}, orient='index', columns=[f'Std T={T}'])
    std_devs[f"Ratio T={T}"] = std_devs[f"Std T={T}"] / T 
    return std_devs

def supercell_convergence_df(temps, supercell_sizes):
    std_temps = []
    for T in temps:
        filepaths = create_filepaths_supercell(supercell_sizes, T)
        data, temperatures = read_supercell_data(filepaths)
        std_devs = get_std_devs(data, supercell_sizes, T)
        std_temps.append(std_devs)
  
    df = pd.concat(std_temps, axis=1)
    df["N_atoms"] = [size**3 * 4 for size in supercell_sizes]
    df.index.name = "Supercell Size"
    return df


def fit_power_law(x, a, b):
    return a * x ** b

def fit_fluctuation_power_law(df, temps, verbose=True, confidence_levels=[1.96]):
    power_law = lmfit.Model(fit_power_law)
    params = power_law.make_params(a=1, b=-1/2)
    fitted_params = {T: None for T in temps}
    fitted_uncertainties = {T: [] for T in temps}
    for T in temps:
        if T == 293: 
            x_data = df["N_atoms"][:-4]
            y_data = df[f"Ratio T={T}"][:-4]
        else:
            x_data = df["N_atoms"]
            y_data = df[f"Ratio T={T}"]
        result = power_law.fit(y_data, params, x=x_data, nan_policy='omit')

        # compute 95% confidence intervals
        fitted_params[T] = result.best_values
        fitted_uncertainties[T].append(result.params['a'].stderr)
        fitted_uncertainties[T].append(result.params['b'].stderr)

        if verbose:
            print(f"\n===== Fit results for T={T} K=====")
            # print(result.fit_report())
            print(f"\nPower law exponent: {result.best_values['b']:.4f} ± {result.params['b'].stderr:.4f}")
            sigma_levels = confidence_levels
            ci = conf_interval(result, result, sigmas=sigma_levels)
            print("\n Confidence Report:")
            report_ci(ci)
            print(f"====================================\n")
    return fitted_params, fitted_uncertainties


# ============================
#        TIMESTEP AND SUPERCELL CONVERGENCE
# ============================

def plot_timestep_and_supercell_convergence(timesteps, filepaths, temps, supercell_sizes, total_energies,
                                            timestep_threshold=1e-5, supercell_threshold=0.03, n_atoms=32, verbose=True,
                                            markers=['o', 'd', 's', 'v', 'D', '+'],
                                            marker_colors=['#005f99', '#cc3300', "#157c3e", '#a64d79', '#ff9933', '#000000'],
                                            line_colors=['#007acc', '#ff5733', '#2ecc71', '#9b59b6', "#ffc972", "#5E5E5E"]):
    std_devs = {t: np.std(total_energies[t])/n_atoms for t in timesteps if t in total_energies}
    df = pd.DataFrame.from_dict(std_devs, orient='index', columns=['Standard Deviation'])
    df["Converged"] = df['Standard Deviation'] < timestep_threshold

    fig, axs = plt.subplots(2,1, figsize=(6, 10))
    fs = 25
    ax = axs[0]
    converged = df[df['Converged']]
    not_converged = df[~df['Converged']]
    if verbose: 
        print("Converged timesteps:")
        print(converged)
        print("\nNot converged timesteps:")
        print(not_converged)
    ax.loglog(converged.index, converged['Standard Deviation'], 'o', color="#8b1616", markersize=10, markerfacecolor="#df3113", markeredgewidth=2, zorder=2)
    ax.loglog(not_converged.index, not_converged['Standard Deviation'], 'o', color="#006ead", markersize=10, markerfacecolor='none', markeredgewidth=2, zorder=1)
    ax.loglog(df.index, df['Standard Deviation'], '-', color='black', alpha=0.5, zorder=1)
    ax.set_xlabel(r"$\boldsymbol{\Delta t}$ $\textbf{[ps]}$", fontsize=fs)
    ax.set_ylabel(r"$\boldsymbol{\sigma_E}$ $\textbf{[eV/atom]}$", fontsize=fs)
    ax.axhline(y=1e-5, color="#000000", linestyle='--', label=rf'Convergence threshold', linewidth=3, zorder=1)
    ax.text(0.83, 0.15, r'$\boldsymbol{(a)}$', transform=ax.transAxes, fontsize=32, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.grid(True)
    ax.legend()

    ax = axs[1]

    df = supercell_convergence_df(temps, supercell_sizes)
    fitted_params, fitted_uncertainties = fit_fluctuation_power_law(df, temps, True, confidence_levels=[1.96])
    for T, marker, marker_color, line_color in zip(temps, markers, marker_colors, line_colors):
        ax.loglog(df["N_atoms"], df[f"Ratio T={T}"], marker=marker, label=rf'${T}$ K', markersize=12, markeredgecolor=marker_color,
                markerfacecolor='none', markeredgewidth=2.5, color=line_color, linestyle='None', zorder=3)

    x_plot = np.linspace(df["N_atoms"].min(), df["N_atoms"].max(), 1000)
    for T, params, uncertainties, color in zip(temps, fitted_params.values(), fitted_uncertainties.values(), line_colors):
        a, b = params['a'], params['b']
        a_err, b_err = uncertainties[0], uncertainties[1]
        
        # b_ci is read from the 95% confidence interval of the fit
        # since we only care about the power law exponent, we consider a as a constant
        b_ci = 0.04 
        y_plot = a * x_plot ** b
        y1 = a * x_plot ** (b - b_ci)
        y2 = a * x_plot ** (b + b_ci)
        ax.fill_between(x_plot, y1, y2, alpha=0.15, color=color, zorder=1)
        ax.loglog(x_plot, y_plot, label=rf'$\propto \boldsymbol N^{{\mathbf{{{b:.2f} \pm {b_ci:.2f}}}}}$',color=color, linewidth=2, zorder=2)
    ax.axhline(y=supercell_threshold, color='black', linestyle='--', linewidth=3, zorder=1)
    ax.set_xlabel(r'$\boldsymbol{N}$', fontsize=fs)
    ax.set_ylabel(r'$\boldsymbol{\sigma_T/T}$', fontsize=fs)

    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=True)
    # plt.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.35), frameon=True)
    ax.legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.05, -0.07))
    ax.grid(which='both', linestyle='-', linewidth=1, alpha=0.2)
    ax.text(0.83, 0.95, r'$\boldsymbol{(b)}$', transform=ax.transAxes, fontsize=32, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    
    plt.tight_layout()

    os.makedirs("Figures", exist_ok=True)
    plt.savefig('Figures/convergence_analysis.pdf', bbox_inches='tight')

    plt.show()


# ============================
#         VELOCITY EQUILIBRATION
# ============================
def load_velocity_distributions(filename):
    df = pd.read_csv(filename, sep=' ', skiprows=9, names=['id', 'type', 'vx', 'vy', 'vz'])
    df['v'] = np.sqrt(df['vx']**2 + df['vy']**2 + df['vz']**2)
    df_I = df[df["type"] == 1].copy()
    df_Ag = df[df["type"] == 2].copy() 
    return df_I, df_Ag

def gaussian(x, temp, mass):
    sigma = np.sqrt(k * temp / mass) / 100  # convertir en Å/ps
    return stats.norm.pdf(x, loc=0, scale=sigma)

def maxwell_boltzmann(x, temp, mass):
    scale = np.sqrt(k * temp / mass) / 100  # convertir en Å/ps
    return stats.maxwell.pdf(x, scale=scale)


def plot_velocity_distributions(df_init, df_equil, T, species='Iodine'):
    if species not in ['Iodine', 'Silver']:
        raise ValueError("Species must be either 'Iodine' or 'Silver'.")

    m_I = 126.90447 * atomic_mass 
    m_Ag = 107.8682 * atomic_mass

    if species == 'Iodine':
        s = 0 
        m = m_I
    elif species == 'Silver':
        s = 1
        m = m_Ag

    # x_values for Gaussian and Maxwell-Boltzmann distributions
    x_min = -7.5 
    x_max = 7.5
    x_vals_gauss = np.linspace(x_min, x_max, 1000)
    x_vals_maxwell = np.linspace(0, 9, 1000)

    fig, axs = plt.subplots(4, 2, figsize=(8, 16))
    n_bins = 15
    hist_kws = {'alpha': 0.3, 'linewidth': 1.5, 'element': 'step', 'stat': 'density', 'bins': n_bins}

    # Theoretical distributions
    for i in range(3): 
        ax = axs[i, 1]
        ax.set_ylim(0, 0.33)
        ax.plot(x_vals_gauss, gaussian(x_vals_gauss, 293, m), '-', color="#0012b6", linewidth=3, label=r'Théorie 293 K')
        ax.plot(x_vals_gauss, gaussian(x_vals_gauss, 800, m), '-', color="#b90000", linewidth=3, label=r'Théorie 800 K')

        ax = axs[i, 0]
        ax.set_ylim(0, 0.33)

    # Maxwell-Boltzmann distribution
    ax = axs[3, 1]
    ax.plot(x_vals_maxwell, maxwell_boltzmann(x_vals_maxwell, 293, m), '-', color="#0012b6", linewidth=3, label=r'Théorie 293 K')
    ax.plot(x_vals_maxwell, maxwell_boltzmann(x_vals_maxwell, 800, m), '-', color="#b90000", linewidth=3, label=r'Théorie 800 K')

    # Plot v_x
    ax = axs[0, 0]
    ax.set_title(r'$\boldsymbol{t=0}$ $\textbf{ps}$', fontsize=32)
    sns.histplot(df_init[293][s]['vx'], ax=ax, label=r'293 K', **hist_kws)
    sns.histplot(df_init[800][s]['vx'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'$\boldsymbol{P(v_x)}$', fontsize=25)
    ax.set_xlabel(r'$\boldsymbol{v_x}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.05, 0.95, r'$\boldsymbol{(a)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(x_min, x_max)

    ax = axs[0, 1]
    ax.set_title(r'$\boldsymbol{t=5}$ $\textbf{ps}$', fontsize=32)
    sns.histplot(df_equil[293][s]['vx'], ax=ax, label=r'293 K', **hist_kws)
    sns.histplot(df_equil[800][s]['vx'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_xlabel(r'$\boldsymbol{v_x}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.8, 0.95, r'$\boldsymbol{(b)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(x_min, x_max)
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0)

    # Plot v_y
    ax = axs[1, 0]
    sns.histplot(df_init[293][s]['vy'], ax=ax, label=r'293 K', **hist_kws)
    sns.histplot(df_init[800][s]['vy'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'$\boldsymbol{P(v_y)}$', fontsize=28)
    ax.set_xlabel(r'$\boldsymbol{v_y}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.05, 0.95, r'$\boldsymbol{(c)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(x_min, x_max)

    ax = axs[1, 1]
    sns.histplot(df_equil[293][s]['vy'], ax=ax, label=r'293 K', **hist_kws)
    sns.histplot(df_equil[800][s]['vy'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_xlabel(r'$\boldsymbol{v_y}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.8, 0.95, r'$\boldsymbol{(d)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(x_min, x_max)
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0)

    # Plot v_z
    ax = axs[2, 0]
    sns.histplot(df_init[293][s]['vz'], ax=ax, label=r'293 K', **hist_kws) 
    sns.histplot(df_init[800][s]['vz'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'$\boldsymbol{P(v_z)}$', fontsize=25)
    ax.set_xlabel(r'$\boldsymbol{v_z}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.05, 0.95, r'$\boldsymbol{(e)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(x_min, x_max)

    ax = axs[2, 1]
    sns.histplot(df_equil[293][s]['vz'], ax=ax, label=r'293 K', **hist_kws)
    sns.histplot(df_equil[800][s]['vz'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_xlabel(r'$\boldsymbol{v_z}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.8, 0.95, r'$\boldsymbol{(f)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(x_min, x_max)
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0)

    # plot |v|
    ax = axs[3, 0]
    sns.histplot(df_init[293][s]['v'], ax=ax, label=r'293 K', **hist_kws)
    sns.histplot(df_init[800][s]['v'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'$\boldsymbol{P(\lVert \vec v\rVert)}$', fontsize=25)
    ax.set_xlabel(r'$\boldsymbol{\lVert \vec v\rVert}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.05, 0.95, r'$\boldsymbol{(g)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 0.75)

    ax = axs[3, 1]
    sns.histplot(df_equil[293][s]['v'], ax=ax, label=r'293 K', **hist_kws)
    sns.histplot(df_equil[800][s]['v'], ax=ax, label=r'800 K', **hist_kws)
    ax.set_ylabel(r'')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_xlabel(r'$\boldsymbol{\lVert \vec v\rVert}$ [$\textbf{\AA/ps}$]')
    ax.grid()
    ax.text(0.8, 0.95, r'$\boldsymbol{(h)}$',  transform=ax.transAxes, fontsize=25, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 0.75)
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0)

    plt.tight_layout()

    legend_elements = [
        Patch(facecolor='tab:blue', label='293 K', alpha=0.5, edgecolor='tab:blue'),
        Patch(facecolor='tab:orange', label='800 K', alpha=0.5, edgecolor='tab:orange'),
        Line2D([0], [0], color='#0012b6', lw=4, linestyle='-', label='Theoretical 293 K'),
        Line2D([0], [0], color='#b90000', lw=4, linestyle='-', label='Theoretical 800 K')
    ]

    # Set xlabel fontsize to 18 for all axes
    for ax in axs.flat:
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.tick_params(axis='both', which='minor', labelsize=18)
        ax.xaxis.label.set_size(25)
        ax.yaxis.label.set_size(25)

    fig.legend(handles=legend_elements, loc='lower center', ncol=2, fontsize=22, frameon=True, bbox_to_anchor=(0.5, -0.08))

    plt.savefig(f"Figures/velocity_distributions_{species}.pdf", bbox_inches='tight')

    plt.show()


# ============================
#         VAF PLOTS
# ============================
def read_VAF(filepath):
    """Reads the VAF data from a file."""
    df = pd.read_csv(filepath, sep='\s+', header=None, skiprows=1)
    df.columns = ['t', 'VAF(t)']
    return df

def plot_VAF(filepaths, T, 
             x_min=0.09, x_max=0.65,
             markers=['o', 'd', 's', 'v', 'D', '+'], 
             marker_colors=['#005f99', '#cc3300', "#157c3e", '#ffcc00', '#9933cc'],
             line_colors=['#007acc', '#ff5733', '#2ecc71', '#ffcc00', '#9933cc']):
    
    list_df_vaf = [read_VAF(filepath) for filepath in filepaths]
    
    fig, ax = plt.subplots(figsize=(7, 5.5))

    # Plot on main axes
    for df_vaf, temp, m, m_color, l_color in zip(list_df_vaf, T, markers, marker_colors, line_colors):
        ax.plot(df_vaf['t'], df_vaf['VAF(t)']/df_vaf['VAF(t)'].iloc[0], 
            label=f'{temp} K', linewidth=2, marker=m, markersize=10,
            markeredgecolor=m_color, markerfacecolor='none', 
            markeredgewidth=2, color=l_color, zorder=3)

    ax.set_xlabel(r'$\boldsymbol{t}$ $\textbf{[ps]}$', fontsize=30)
    ax.set_ylabel(r'$\textbf{NVAF}$', fontsize=30)
    ax.set_xlim(0, 1.5)
    ax.grid()
    ax.legend(ncol=2)

    # Create inset axes
    axins = inset_axes(ax, width="60%", height="35%", loc='center right', borderpad=0.5)
    axins.spines['top'].set_linewidth(0)
    axins.spines['right'].set_linewidth(0)
    axins.spines['bottom'].set_linewidth(0.5)
    axins.spines['left'].set_linewidth(0.5)
    axins.tick_params(axis='both', labelsize=15, width=0.5, length=4)

    y_min, y_max = -0.09, 0.05

    for df_vaf, temp, m, m_color, l_color in zip(list_df_vaf, T, markers, marker_colors, line_colors):
        axins.plot(df_vaf['t'], df_vaf['VAF(t)']/df_vaf['VAF(t)'].iloc[0],
                linewidth=1.5, marker=m, markersize=8,
                markeredgecolor=m_color, markerfacecolor='none',
                markeredgewidth=1.5, color=l_color)
    axins.set_xlim(x_min, x_max)
    axins.set_ylim(y_min, y_max)
    axins.grid(True, alpha=0.2)

    # Optional: Mark the inset area on the main plot
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.3")

    plt.tight_layout()

    os.makedirs("Figures", exist_ok=True)
    plt.savefig("Figures/vaf_plot.pdf", bbox_inches='tight')

    plt.show()

    return list_df_vaf


def produce_csv_for_vaf(json_file):
    # Lecture du fichier JSON
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Extraction des données
    velocities_Ag = data['velocities_Ag']
    velocities_I = data['velocities_I']
    times = data['times']

    # For Ag
    all_ag_vels = []
    for time_idx, time_val in enumerate(times):
        for atom_idx in range(len(velocities_Ag[0])):
            all_ag_vels.append({
                'time': time_val,
                'time_step': time_idx,
                'atom_idx': atom_idx,
                'vx [m/s]': velocities_Ag[time_idx][atom_idx][0] * 100, # Convert A/ps to m/s
                'vy [m/s]': velocities_Ag[time_idx][atom_idx][1] * 100,
                'vz [m/s]': velocities_Ag[time_idx][atom_idx][2] * 100,
                'v_magnitude [m/s]': np.sqrt(sum([(v*100)**2 for v in velocities_Ag[time_idx][atom_idx]]))
            })

    # For I
    all_i_vels = []
    for time_idx, time_val in enumerate(times):
        for atom_idx in range(len(velocities_I[0])):
            all_i_vels.append({
                'time': time_val,
                'time_step': time_idx,
                'atom_idx': atom_idx,
                'vx [m/s]': velocities_I[time_idx][atom_idx][0] * 100,
                'vy [m/s]': velocities_I[time_idx][atom_idx][1] * 100,
                'vz [m/s]': velocities_I[time_idx][atom_idx][2] * 100,
                'v_magnitude [m/s]': np.sqrt(sum([(v*100)**2 for v in velocities_I[time_idx][atom_idx]]))
            })

    df_Ag = pd.DataFrame(all_ag_vels)
    df_I = pd.DataFrame(all_i_vels)

    m_Ag = 107.8682 * atomic_mass # in kg
    m_I = 126.90447 * atomic_mass 

    df_Ag['kinetic_energy [J]'] = 0.5 * m_Ag * (df_Ag['vx [m/s]']**2 + df_Ag['vy [m/s]']**2 + df_Ag['vz [m/s]']**2) # in kg m^2/s^2 = J
    df_I['kinetic_energy [J]'] = 0.5 * m_I * (df_I['vx [m/s]']**2 + df_I['vy [m/s]']**2 + df_I['vz [m/s]']**2)


    # Save the DataFrames to CSV files
    df_Ag.to_csv('P1/VAF/velocities_Ag_6_800_0.001.csv', index=False)
    df_I.to_csv('P1/VAF/velocities_I_6_800_0.001.csv', index=False)

    return df_Ag, df_I


# ============================
#         RDF PLOTS
# ============================
def load_rdf_data(filepath):
    return pd.read_csv(filepath, sep='\s+', skiprows=1, names=["r", "RDF_I", "IRDF_I", "RDF_Ag", "IRDF_Ag"])

def plot_rdf(list_df_rdf_293, list_df_rdf_10, n_bins, line_colors=['#007acc', '#ff5733', '#2ecc71', '#ffcc00', '#9933cc']):
    fig, axs = plt.subplots(1, 2, figsize=(9, 5)) 
    for ax in axs:
        ax.tick_params(axis='both', which='major', labelsize=25)
        ax.tick_params(axis='both', which='minor', labelsize=25)

    ax = axs[0]
    if len(n_bins) == 1:
        df_rdf = list_df_rdf_293[0]    
        ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I', linewidth=3, color=line_colors[0])
        ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag', linewidth=3, color=line_colors[1])
        max_rdf = df_rdf['RDF_I'].max()
        max_r = df_rdf['r'].max()
        ax.set_xlim(0, max_r)
        ax.set_yticks([0,1,2,3,4])
        ax.set_xticks([0,5,10,15])
        ax.set_ylim(0, 1.05 * max_rdf)
    else:
        for df_rdf, n in zip(list_df_rdf_293, n_bins):
            ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I {n} bins', linewidth=3)
            ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag {n} bins', linewidth=3)
        max_rdf = np.max([df_rdf['RDF_I'].max() for df_rdf in list_df_rdf_293])
        max_r = np.max([df_rdf['r'].max() for df_rdf in list_df_rdf_293])
        ax.set_xlim(0, max_r)
        ax.set_ylim(0, 1.05 * max_rdf)
    ax.text(0.03, 0.95, r'$\boldsymbol{(a)}$', transform=ax.transAxes, fontsize=30, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlabel(r'$\boldsymbol{r}$ $[\textbf{\AA}]$', fontsize=30)
    ax.set_ylabel(r'$\boldsymbol{g(r)}$', fontsize=30)
    ax.legend(title=r"$\boldsymbol{T=293}$ $\textbf{K}$", title_fontsize=25)
    ax.grid()

    ax = axs[1]
    if len(n_bins) == 1:
        df_rdf = list_df_rdf_10[0]    
        ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I', linewidth=3, color=line_colors[0])
        ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag', linewidth=3, color=line_colors[1])
        max_rdf = df_rdf['RDF_I'].max()
        max_r = df_rdf['r'].max()
        ax.set_xlim(0, max_r)
        ax.set_ylim(0, 1.05 * max_rdf)
        ax.set_xticks([0,5,10,15])
    else:
        for df_rdf, n in zip(list_df_rdf_10, n_bins):
            ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I {n} bins', linewidth=3)
            ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag {n} bins', linewidth=3)
        max_rdf = np.max([df_rdf['RDF_I'].max() for df_rdf in list_df_rdf_10])
        max_r = np.max([df_rdf['r'].max() for df_rdf in list_df_rdf_10])
        ax.set_xlim(0, max_r)
        ax.set_ylim(0, 1.05 * max_rdf)
    ax.set_xlabel(r'$\boldsymbol{r}$ $[\textbf{\AA}]$', fontsize=30)
    ax.set_ylabel(r'')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.text(0.03, 0.95, r'$\boldsymbol{(b)}$', transform=ax.transAxes, fontsize=30, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.legend(title=r"$\boldsymbol{T=10}$ $\textbf{K}$", title_fontsize=25, loc='upper right')
    ax.grid()

    plt.tight_layout()

    os.makedirs("Figures", exist_ok=True)
    plt.savefig("Figures/rdf_plot.pdf", bbox_inches='tight')

    plt.show()


# ===============================
#         MSD PLOTS AND FITS
# ===============================
def plot_msd_AgI(list_df_msd, T, t_sim, dt, 
                 line_colors=['#007acc', '#ff5733', '#2ecc71', '#9b59b6', "#ffc972", "#5E5E5E", 
                              "#8b1616", "#df3113", "#0012b6", "#b90000"]):
    x_ticks = np.arange(0, t_sim + dt, dt)

    max_Ag = max([df["MSD_Ag"].max() for df in list_df_msd])
    max_I = max([df["MSD_I"].max() for df in list_df_msd])

    fig, ax = plt.subplots(1, 2, figsize=(10, 6))
    for a in ax:
        a.tick_params(axis='both', which='major', labelsize=25)
        a.tick_params(axis='both', which='minor', labelsize=25)

    for df_msd, T, l_color in zip(list_df_msd, T, line_colors):
        ax[0].plot(df_msd["t (ps)"]-df_msd["t (ps)"].iloc[0], df_msd["MSD_I"], color=l_color, linewidth=4)
        ax[1].plot(df_msd["t (ps)"]-df_msd["t (ps)"].iloc[0], df_msd["MSD_Ag"], label=f"{T} K", color=l_color, linewidth=4)
    ax[0].set_title(r'$\textbf{Iodine}$', fontsize=34)
    ax[0].set_xlabel(r'$\boldsymbol{t}$ $\textbf{[ps]}$', fontsize=32)
    ax[0].set_ylabel(r'$\textbf{MSD}$ [$\textbf{\AA}^2$]', fontsize=32)
    ax[0].grid()
    ax[0].set_xlim(0, t_sim)
    ax[0].set_ylim(0, max_I * 1.03)
    ax[0].set_xticks(x_ticks)
    ax[0].text(0.05, 0.95, r'$\boldsymbol{(a)}$', transform=ax[0].transAxes, fontsize=34, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))

    ax[1].grid()
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].set_xlabel(r'$\boldsymbol{t}$ $\textbf{[ps]}$', fontsize=32)
    ax[1].set_title(r'$\textbf{Silver}$', fontsize=34)
    ax[1].set_xticks(x_ticks)
    ax[1].set_xlim(0, t_sim)
    ax[1].set_ylim(0, max_Ag * 1.03)
    ax[1].text(0.05, 0.95, r'$\boldsymbol{(b)}$', transform=ax[1].transAxes, fontsize=34, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))

    handles, labels = ax[1].get_legend_handles_labels()
    # ax[1].legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.15, 0.5), frameon=True, fontsize=30)
    
    fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=25, frameon=True, bbox_to_anchor=(0.51, -0.22))

    plt.tight_layout()
    plt.savefig("Figures/msd_AgI_plot.pdf", bbox_inches='tight')
    plt.show()


def fit_msd(list_df_msd, T, model, t_min=0, confidence_levels=[1.96], verbose=True):
    slopes = {t: None for t in T}
    uncertainties = {t: None for t in T}   
    results = {t: None for t in T} 
    for df_msd, temp in zip(list_df_msd, T):
        mask = (df_msd["t (ps)"] - df_msd["t (ps)"].iloc[0]) >= t_min
        t = df_msd["t (ps)"][mask] - df_msd["t (ps)"].iloc[0] 
        msd_Ag = df_msd["MSD_Ag"][mask]
        params = model.make_params(a=2, b=1)
        result = model.fit(msd_Ag, params, x=t)
        a = result.best_values['a']
        a_err = result.params['a'].stderr
        slopes[temp] = a
        uncertainties[temp] = a_err
        results[temp] = result
        if verbose:
            print(f"==== Results for T={temp} K ====")
            print(f"slope: {a:.4f} ± {a_err:.4f} (units: ps^(-1)*A^2)")                  
            sigma_levels = confidence_levels
            ci = conf_interval(result, result, sigmas=sigma_levels)
            print("\n Confidence Report:")
            report_ci(ci)
            print("\n")
    df_results = pd.DataFrame({
        'Temperature (K)': T,
        'Slope (ps^(-1)*A^2)': [slopes[t] for t in T],
        'Uncertainty (ps^(-1)*A^2)': [uncertainties[t] for t in T]
    })
    return df_results, results

def linear_function(x, a, b):
    return a * x + b
