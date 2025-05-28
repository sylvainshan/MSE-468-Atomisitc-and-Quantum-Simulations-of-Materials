"""
Reads the output json file, and make an XCRYSDEN-readable axsf file that can be visualised.
By Leonid Kahle (Spring 2018)
"""

import sys, json

def create_axsf(cell, positions_data: dict):
    n_images = len(positions_data[list(positions_data.keys())[0]])
    n_atoms = len(positions_data[list(positions_data.keys())[0]][0])*len(positions_data)
    text = 'ANIMSTEPS %s \n' % ( n_images )
    text += """CRYSTAL
PRIMVEC
%12.8f   %12.8f   %12.8f
%12.8f   %12.8f   %12.8f
%12.8f   %12.8f   %12.8f
""" % ( cell[0][0],cell[0][1],cell[0][2],cell[1][0],cell[1][1],cell[1][2],cell[2][0],cell[2][1],cell[2][2] )
    for i_image in range(n_images):
        text += 'PRIMCOORD %s\n' % (i_image+1)
        text += '%s 1\n' % ( n_atoms )
        for element, positions in positions_data.items():
            for this_atom_pos in positions[i_image]:
                text += '%s  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n' % (element, this_atom_pos[0],this_atom_pos[1],this_atom_pos[2],0,0,0 )

    return text

def read_positions_cell(json_file):
    with open(json_file,'r') as f:
        all_dict = json.load(f)

    return {"Ag": all_dict['positions_Ag'],
            "I": all_dict['positions_I']}, all_dict['cell']

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('json_file', help='The json output, produces with parser.py'
        ' from a gulp trajectory')
    parser.add_argument('-o', '--output-axsf', help='Filename to write animated XSF file', default='xcrysden.axsf')

    parsed_args = parser.parse_args()

    positions, cell = read_positions_cell(parsed_args.json_file)

    axsf_text = create_axsf(cell, positions)

    with open(parsed_args.output_axsf,'w') as f:
        f.write(axsf_text)
