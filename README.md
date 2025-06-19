# Materials Simulations with DFT and Molecular Dynamics

This repository contains a series of simulation projects (`lab1`, `lab2`, `lab3`, `lab4`) focused on materials modeling using **Density Functional Theory (DFT)** and **Molecular Dynamics (MD)**. These simulations are performed with industry-standard tools such as **Quantum ESPRESSO** and **LAMMPS**.

Each lab directory corresponds to a separate project, covering different physical systems or simulation objectives.

## 🧪 Contents

- **Lab 1 – FCC Silver Properties with LAMMPS**  
  - Structural relaxation with Lennard-Jones and EAM potentials  
  - Bulk modulus estimation via energy–volume fitting  
  - Vacancy formation energy (with and without relaxation)  
  - Surface energy calculation 
  - Modeling of the $\rm L1_0$ CuAu binary alloy  

- **Lab 2 – First-Principles Study of Bulk CaO**  
  DFT simulations of bulk calcium oxide. This lab includes:
  - Convergence tests for energy cutoff and $\mathbf{k}$-point grid  
  - Determination of equilibrium lattice parameter  
  - Calculation of bulk modulus via Birch-Murnaghan fitting  
  - Extraction of elastic constants ($C_{11}$, $C_{12}$, $C_{44}$) from strain–energy relations  

- **Lab 3 – Electronic Properties from First-Principles**  
  - Band structure and density of states of FCC MgO  
  - van der Waals effects on structural and electronic properties of bulk and monolayer $\text{MoS}_2$  
  - Magnetic phase stability of FCC and HCP cobalt, including spin-polarized density of states analysis  

- **Lab 4 – Classical Molecular Dynamics of Silver Iodide**
  - Convergence tests (timestep, supercell size, velocity autocorrelation function decay)  
  - Structural analysis via radial distribution function  
  - Diffusion behavior via mean square displacement 
  - Temperature-dependent diffusion coefficient and Arrhenius fit 

## 📂 Reports

All written reports summarizing the methods, results, and conclusions for each lab are stored in the `reports/` directory. These are provided as PDF files and complement the simulations and notebook analyses.

## 🧰 Tools & Frameworks

- **Quantum ESPRESSO** – Used for electronic structure calculations based on DFT.
- **LAMMPS** – Used for classical molecular dynamics simulations.
- **Jupyter Notebooks** – All analysis and visualization are performed in interactive Jupyter notebooks.

## 📈 Analysis

Each lab includes Jupyter notebooks for:
- Pre-processing and input preparation
- Post-processing and data analysis
- Visualization of physical properties (e.g., energy, temperature, radial distribution functions, etc.)
