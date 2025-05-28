import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import os
import lmfit
from lmfit import report_ci, conf_interval, Model
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


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
    
    plt.figure(figsize=(6, 5))
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
#         VAF PLOTS
# ============================
def read_VAF(filepath):
    """Reads the VAF data from a file."""
    df = pd.read_csv(filepath, sep='\s+', header=None, skiprows=1)
    df.columns = ['t', 'VAF(t)']
    return df

def plot_VAF(filepaths, T, 
             markers=['o', 'd', 's', 'v', 'D', '+'], 
             marker_colors=['#005f99', '#cc3300', "#157c3e", '#ffcc00', '#9933cc'],
             line_colors=['#007acc', '#ff5733', '#2ecc71', '#ffcc00', '#9933cc']):
    
    list_df_vaf = [read_VAF(filepath) for filepath in filepaths]
    
    fig, ax = plt.subplots(figsize=(6, 5))
    # Plot on main axes
    for df_vaf, temp, m, m_color, l_color in zip(list_df_vaf, T, markers, marker_colors, line_colors):
        ax.plot(df_vaf['t'], df_vaf['VAF(t)']/df_vaf['VAF(t)'].iloc[0], 
            label=f'{temp} K', linewidth=2, marker=m, markersize=10,
            markeredgecolor=m_color, markerfacecolor='none', 
            markeredgewidth=2, color=l_color, zorder=3)

    ax.set_xlabel(r'$\boldsymbol{t}$ [ps]')
    ax.set_ylabel(r'$\textbf{VAF}\boldsymbol{(t)}$')
    ax.set_xlim(0, 1.5)
    ax.grid()
    ax.legend(ncol=2)

    # Create inset axes
    axins = inset_axes(ax, width="60%", height="40%", loc='center right', borderpad=0.5)
    axins.spines['top'].set_linewidth(0)
    axins.spines['right'].set_linewidth(0)
    axins.spines['bottom'].set_linewidth(0.5)
    axins.spines['left'].set_linewidth(0.5)
    axins.tick_params(axis='both', labelsize=15, width=0.5, length=4)

    x_min, x_max = 0.09, 0.65
    y_min, y_max = -0.05, 0.05

    for df_vaf, temp, m, m_color, l_color in zip(list_df_vaf, T, markers, marker_colors, line_colors):
        axins.plot(df_vaf['t'], df_vaf['VAF(t)']/df_vaf['VAF(t)'].iloc[0],
                linewidth=1.5, marker=m, markersize=8,
                markeredgecolor=m_color, markerfacecolor='none',
                markeredgewidth=1.5, color=l_color)
    axins.set_xlim(x_min, x_max)
    axins.set_ylim(y_min, y_max)
    axins.grid(True, alpha=0.2)

    # Optional: Mark the inset area on the main plot
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    plt.tight_layout()

    os.makedirs("Figures", exist_ok=True)
    plt.savefig("Figures/vaf_plot.pdf", bbox_inches='tight')

    plt.show()

    return list_df_vaf


# ============================
#         RDF PLOTS
# ============================
def load_rdf_data(filepath):
    return pd.read_csv(filepath, sep='\s+', skiprows=1, names=["r", "RDF_I", "IRDF_I", "RDF_Ag", "IRDF_Ag"])

def plot_rdf(list_df_rdf_293, list_df_rdf_10, n_bins, line_colors=['#007acc', '#ff5733', '#2ecc71', '#ffcc00', '#9933cc']):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.5)) 
    ax = axs[0]
    if len(n_bins) == 1:
        df_rdf = list_df_rdf_293[0]    
        ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I', linewidth=2, color=line_colors[0])
        ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag', linewidth=2, color=line_colors[1])
        max_rdf = df_rdf['RDF_I'].max()
        max_r = df_rdf['r'].max()
        ax.set_xlim(0, max_r)
        ax.set_ylim(0, 1.05 * max_rdf)
    else:
        for df_rdf, n in zip(list_df_rdf_293, n_bins):
            ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I {n} bins', linewidth=2)
            ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag {n} bins', linewidth=2)
        max_rdf = np.max([df_rdf['RDF_I'].max() for df_rdf in list_df_rdf_293])
        max_r = np.max([df_rdf['r'].max() for df_rdf in list_df_rdf_293])
        ax.set_xlim(0, max_r)
        ax.set_ylim(0, 1.05 * max_rdf)
    ax.text(0.03, 0.95, r'$\boldsymbol{(a)}$', transform=ax.transAxes, fontsize=30, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.set_xlabel(r'$\boldsymbol{r}$ $[\textbf{\AA}]$')
    ax.set_ylabel(r'$\boldsymbol{g(r)}$')
    ax.legend(title=r"$T=293$ K", title_fontsize=22)
    ax.grid()

    ax = axs[1]
    if len(n_bins) == 1:
        df_rdf = list_df_rdf_10[0]    
        ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I', linewidth=2, color=line_colors[0])
        ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag', linewidth=2, color=line_colors[1])
        max_rdf = df_rdf['RDF_I'].max()
        max_r = df_rdf['r'].max()
        ax.set_xlim(0, max_r)
        ax.set_ylim(0, 1.05 * max_rdf)
    else:
        for df_rdf, n in zip(list_df_rdf_10, n_bins):
            ax.plot(df_rdf['r'], df_rdf['RDF_I'], label=f'I {n} bins', linewidth=2)
            ax.plot(df_rdf['r'], df_rdf['RDF_Ag'], label=f'Ag {n} bins', linewidth=2)
        max_rdf = np.max([df_rdf['RDF_I'].max() for df_rdf in list_df_rdf_10])
        max_r = np.max([df_rdf['r'].max() for df_rdf in list_df_rdf_10])
        ax.set_xlim(0, max_r)
        ax.set_ylim(0, 1.05 * max_rdf)
    ax.set_xlabel(r'$\boldsymbol{r}$ $[\textbf{\AA}]$')
    ax.set_ylabel(r'')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.text(0.03, 0.95, r'$\boldsymbol{(b)}$', transform=ax.transAxes, fontsize=30, verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='white', alpha=1, edgecolor='None'))
    ax.legend(title=r"$T=10$ K", title_fontsize=22)
    ax.grid()

    plt.tight_layout()

    os.makedirs("Figures", exist_ok=True)
    plt.savefig("Figures/rdf_plot.pdf", bbox_inches='tight')

    plt.show()


# ===============================
#         MSD PLOTS AND FITS
# ===============================
def plot_msd_AgI(list_df_msd, T, t_sim, dt, line_colors=['#007acc', '#ff5733', '#2ecc71', '#9b59b6', "#ffc972", "#5E5E5E"]):
    x_ticks = np.arange(0, t_sim + dt, dt)

    max_Ag = max([df["MSD_Ag"].max() for df in list_df_msd])
    max_I = max([df["MSD_I"].max() for df in list_df_msd])

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for df_msd, T, l_color in zip(list_df_msd, T, line_colors):
        ax[0].plot(df_msd["t (ps)"]-df_msd["t (ps)"].iloc[0], df_msd["MSD_I"], color=l_color, linewidth=2)
        ax[1].plot(df_msd["t (ps)"]-df_msd["t (ps)"].iloc[0], df_msd["MSD_Ag"], label=f"{T} K", color=l_color, linewidth=2)
    ax[0].set_title(r'$\textbf{Iodine}$')
    ax[0].set_xlabel(r'$\boldsymbol{t}$ $\textbf{[ps]}$')
    ax[0].set_ylabel(r'$\textbf{MSD}$ [$\textbf{\AA}^2$]')
    ax[0].grid()
    ax[0].set_xlim(0, t_sim)
    ax[0].set_ylim(0, max_I * 1.03)
    ax[0].set_xticks(x_ticks)

    ax[1].grid()
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].set_xlabel(r'$\boldsymbol{t}$ $\textbf{[ps]}$')
    ax[1].set_title(r'$\textbf{Silver}$')
    ax[1].set_xticks(x_ticks)
    ax[1].set_xlim(0, t_sim)
    ax[1].set_ylim(0, max_Ag * 1.03)

    handles, labels = ax[1].get_legend_handles_labels()
    ax[1].legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.1, 0.5), frameon=True)

    plt.tight_layout()
    plt.show()


def fit_msd(list_df_msd, T, model, t_min=0, confidence_levels=[1.96], verbose=True):
    slopes = {t: None for t in T}
    uncertainties = {t: None for t in T}    
    for df_msd, temp in zip(list_df_msd, T):
        mask = (df_msd["t (ps)"] - df_msd["t (ps)"].iloc[0]) >= t_min
        t = df_msd["t (ps)"][mask] - df_msd["t (ps)"].iloc[0] 
        msd_Ag = df_msd["MSD_Ag"][mask]
        params = model.make_params(a=1, b=0.0)
        result = model.fit(msd_Ag, params, x=t)
        a = result.best_values['a']
        a_err = result.params['a'].stderr
        slopes[temp] = a
        uncertainties[temp] = a_err
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
    return df_results

def linear_function(x, a, b):
    return a * x + b
