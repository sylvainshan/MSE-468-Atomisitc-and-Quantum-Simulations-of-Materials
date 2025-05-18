import matplotlib 
import matplotlib.pyplot as plt
import numpy as np

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