import os

# Define the directory containing the XYZ files
directory_path = '../'

# Define the header for the Quick input files
header = """HF BASIS=STO-3G cutoff=1.0e-10 denserms=1.0e-6 OPTIMIZE

"""

def convert_xyz_to_quick_input(directory):
    for filename in os.listdir(directory):
        if filename.endswith('.xyz'):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as xyz_file:
                xyz_content = xyz_file.readlines()

            # Remove the first two lines which contain the atom count and comment
            xyz_coords = xyz_content[2:]

            # Create the Quick input file content
            quick_input_content = header + ''.join(xyz_coords)

            # Write the content to a new file with .inp extension
            quick_input_filename = os.path.splitext(filename)[0] + '.in'
            quick_input_file_path = os.path.join(directory, quick_input_filename)
            with open(quick_input_file_path, 'w') as quick_input_file:
                quick_input_file.write(quick_input_content)

# Run the conversion function
convert_xyz_to_quick_input(directory_path)

