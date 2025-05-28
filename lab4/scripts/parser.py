#!/usr/bin/python
"""
Reads the LAMMPS outputs and convert the data to a json file that is used for all
postprocessing steps in LAB 4
By Tushar Thakur (Spring 2024)
"""

import sys, json

def parse_lammps_output(output_fname, positions_fname, velocities_fname, json_fname, e_of_t_fname):
    """
    :param output_fname: The name of the output file produced by LAMMPS
    :param positions_fname: The name of the positions trajectory produced by LAMMPS
    :param velocities_fname: The name of the velocities trajectory produced by LAMMPS
    :param json_fname: The filename I will dump my result to, as a json
    :param e_of_t_fname: The filename I will write extensive properties of the system to
    """

    if not json_fname.endswith('.json'):
        json_fname+='.json'

    with open(positions_fname,'r') as f:
        pos_lines = f.readlines()

    with open(velocities_fname,'r') as f:
        vel_lines = f.readlines()

    with open(output_fname, 'r') as f:
        outlines = f.readlines()

    positions_Ag = []
    velocities_Ag = []
    positions_I = []
    velocities_I = []

    for count, line in enumerate(pos_lines, start=1):
        if "ITEM: NUMBER OF ATOMS" in line:
            nat = int(pos_lines[count])
        if "ITEM: ATOMS id type xu yu zu" in line:
            frame_data_Ag = []
            frame_data_I = []
            for i in range(nat):
                atom_data = pos_lines[count + i].split()
                atom_type = int(atom_data[1])
                x = float(atom_data[2])
                y = float(atom_data[3])
                z = float(atom_data[4])
                if atom_type == 1:
                    frame_data_I.append([x, y, z])
                elif atom_type == 2:
                    frame_data_Ag.append([x, y, z])
            positions_Ag.append(frame_data_Ag)
            positions_I.append(frame_data_I)
    positions_Ag = positions_Ag[1:]
    positions_I = positions_I[1:]

    for count, line in enumerate(vel_lines, start=1):
        if "ITEM: NUMBER OF ATOMS" in line:
            nat = int(vel_lines[count])
        if "ITEM: ATOMS id type vx vy vz" in line:
            frame_data_I = []
            frame_data_Ag = []
            for i in range(nat):
                atom_data = vel_lines[count + i].split()
                atom_type = int(atom_data[1])
                x = float(atom_data[2])
                y = float(atom_data[3])
                z = float(atom_data[4])
                if atom_type == 1:
                    frame_data_I.append([x, y, z])
                elif atom_type == 2:
                    frame_data_Ag.append([x, y, z])
            velocities_Ag.append(frame_data_Ag)
            velocities_I.append(frame_data_I)
    velocities_I = velocities_I[1:]
    velocities_Ag = velocities_Ag[1:]

    for count, line in enumerate(outlines, start=1):
        if "Time           Temp          PotEng         KinEng         Press" in line:
            line_to_start = count+1
            break

    line_to_end = line_to_start + len(positions_Ag)

    times = []
    temperatures = []
    potential_energies = []
    kinetic_energies = []
    pressures = []
    msd_I = []
    msd_Ag = []
    for line in outlines[line_to_start:line_to_end]:
        thermo = line.split()
        times.append(float(thermo[0]))
        temperatures.append(float(thermo[1]))
        potential_energies.append(float(thermo[2]))
        kinetic_energies.append(float(thermo[3]))
        pressures.append(float(thermo[4]))
        msd_I.append(float(thermo[5]))
        msd_Ag.append(float(thermo[6]))

    for count, line in enumerate(outlines, start=1):
        if "Cella          Cellb          Cellc" in line:
            line_to_start = count
            break

    for line in outlines[line_to_start:]:
        abc = line.split()
        try:
            cell = [[float(abc[0]), 0, 0], [0, float(abc[1]), 0], [0, 0, float(abc[2])]]
        except ValueError: 
            break

    all_dict = {}
    all_dict['velocities_Ag'] = velocities_Ag
    all_dict['velocities_I'] = velocities_I
    all_dict['positions_Ag'] = positions_Ag
    all_dict['positions_I'] = positions_I
    all_dict['times'] = times
    all_dict['temperatures'] = temperatures
    all_dict['kinetic_energies'] = kinetic_energies
    all_dict['potential_energies'] = potential_energies
    all_dict['pressures'] = pressures
    all_dict['msd_I'] = msd_I
    all_dict['msd_Ag'] = msd_Ag
    all_dict['cell'] = cell

    with open(json_fname,'w') as f:
        json.dump(all_dict,f)

    with open(e_of_t_fname,'w') as f:
        f.write('# \t t (ps)\t T (K)\tE (eV)\tK (eV)\tP (Bar)\tMSD_I(Å²)\tMSD_Ag(Å²)\n')
        for i in zip(times, temperatures, potential_energies, kinetic_energies, pressures, msd_I, msd_Ag):
            f.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (i[0],i[1],i[2],i[3],i[4],i[5],i[6]) )


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output', help='The output file produced by LAMMPS')
    parser.add_argument('positions', help='The trajectory file containing positions produced by LAMMPS')
    parser.add_argument('velocities', help='The trajectory file containing positions produced by LAMMPS')
    parser.add_argument('json_output', help='The name of the JSON-output')
    parser.add_argument('-e', '--e-of-t', default='E_of_t.dat',
        help='The filename to write energy and pressure as a function of time, default: E_of_t.dat')

    parsed_args = parser.parse_args()
    parse_lammps_output(parsed_args.output, parsed_args.positions, parsed_args.velocities,  parsed_args.json_output, parsed_args.e_of_t)
