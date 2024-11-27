import numpy as np

# Function to read the x,y data from a txt file
def read_xy_file(file_path):
    x_data = []
    y_data = []
    
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line by ',' and strip whitespace
            x, y = line.strip().split()
            x_data.append(float(x))
            y_data.append(float(y))
    
    return x_data, y_data