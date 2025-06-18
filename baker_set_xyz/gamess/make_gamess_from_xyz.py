import os

# Define the directory containing the XYZ files
directory_path = '../'

# Define the header for the GAMESS input files
header = """!   File created by Python script
 $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE MAXIT=200 MULT=1 $END
 $SYSTEM mwords=1000 $END
 $contrl nprint=-5 $end
 $BASIS GBASIS=STO NGAUSS=3 $END
 $SCF DIRSCF=.TRUE. $END
 $STATPT OPTTOL=0.0003 NSTEP=200 method=RFO $END
 $DATA
Title
C1
"""

def convert_xyz_to_gamess_input(directory):
    for filename in os.listdir(directory):
        if filename.endswith('.xyz'):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as xyz_file:
                xyz_content = xyz_file.readlines()

            # Remove the first two lines which contain the atom count and comment
            xyz_coords = xyz_content[2:]

            # Convert XYZ coordinates to GAMESS format
            gamess_coords = []
            for line in xyz_coords:
                parts = line.split()
                element = parts[0]
                x, y, z = parts[1], parts[2], parts[3]
                atomic_number = get_atomic_number(element)
                gamess_coords.append(f" {element} {atomic_number}.0    {x}    {y}    {z}")

            # Create the GAMESS input file content
            gamess_input_content = header + '\n'.join(gamess_coords) + "\n $END\n"

            # Write the content to a new file with .inp extension
            gamess_input_filename = os.path.splitext(filename)[0] + '.inp'
            gamess_input_file_path = os.path.join(directory, gamess_input_filename)
            with open(gamess_input_file_path, 'w') as gamess_input_file:
                gamess_input_file.write(gamess_input_content)

def get_atomic_number(element):
    # Define a dictionary mapping elements to their atomic numbers
    atomic_numbers = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
        'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
        # Add more elements as needed
    }
    return atomic_numbers.get(element, 0)

# Run the conversion function
convert_xyz_to_gamess_input(directory_path)

