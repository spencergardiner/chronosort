import numpy as np
import sys

# Parse the RMSD data from the XPM text file
def parse_rmsd_data(filename):
    # Read the entire file
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Extract the dimensions from the first line
    dimensions_line = lines[0].strip()
    parts = dimensions_line.split()
    if len(parts) >= 2:
        try:
            n_rows = int(parts[0])
            n_cols = int(parts[1])
            print(f"Matrix dimensions: {n_rows}x{n_cols}")
        except ValueError:
            print("Warning: Could not parse matrix dimensions, using fallback")
            n_rows = n_cols = 0  # Will be determined by counting
    else:
        n_rows = n_cols = 0  # Will be determined by counting
    
    # Find where the color definitions end and the matrix data begins
    matrix_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith('"') and not any(c in line for c in ['c #', '/*']):
            matrix_start = i
            break
    
    # Create a color to value mapping from the color definitions
    color_map = {}
    for i in range(1, matrix_start):
        line = lines[i].strip()
        if 'c #' in line and '/*' in line:
            parts = line.split()
            if len(parts) >= 4:
                color_idx = int(parts[0])
                value_str = parts[-1].strip('/*').strip('*/')
                try:
                    color_map[color_idx] = float(value_str)
                except ValueError:
                    print(f"Warning: Could not parse RMSD value: {line}")
    
    # Parse the matrix data
    matrix_data = []
    for i in range(matrix_start, len(lines)):
        line = lines[i].strip()
        if not line.startswith('"'):
            continue
        
        # Extract the row data string inside the quotes
        row_str = line.split('"')[1] if '"' in line else ""
        if not row_str:
            continue
        
        # Convert each character to its corresponding RMSD value
        row_data = [color_map.get(int(c), 0.0) if c.isdigit() else 0.0 for c in row_str]
        matrix_data.append(row_data)
    
    # Check dimensions
    if n_rows > 0 and len(matrix_data) != n_rows:
        print(f"Warning: Expected {n_rows} rows but found {len(matrix_data)}")
    
    if len(matrix_data) > 0 and n_cols > 0 and len(matrix_data[0]) != n_cols:
        print(f"Warning: Expected {n_cols} columns but found {len(matrix_data[0])}")
    
    return np.array(matrix_data)

# Get the path from command line
rmsd_file = sys.argv[1]
output_file = sys.argv[2]
frame_count = int(sys.argv[3])

# Parse RMSD data
rmsd_matrix = parse_rmsd_data(rmsd_file)

# Make sure we have a square matrix
n = rmsd_matrix.shape[0]
if n < frame_count:
    print(f"Error: Expected at least {frame_count} frames, but found {n}")
    sys.exit(1)

# Create mask for visited frames (1 = visited)
mask = np.zeros(n, dtype=bool)

# Start with a random frame
start_frame = np.random.randint(0, n)
ordered_frames = [start_frame]
mask[start_frame] = True

# Build the ordered list
while len(ordered_frames) < frame_count:
    last_frame = ordered_frames[-1]
    
    # Get the row of RMSD values for the last frame
    rmsd_row = rmsd_matrix[last_frame]
    
    # Create a masked version where visited frames have RMSD = infinity
    masked_rmsd = np.where(mask, np.inf, rmsd_row)
    
    # Find the frame with minimum RMSD to the last frame
    next_frame = np.argmin(masked_rmsd)
    
    # Add to our ordered list and mark as visited
    ordered_frames.append(next_frame)
    mask[next_frame] = True

# Output the ordered frame indices
with open(output_file, 'w') as f:
    for idx in ordered_frames:
        f.write(f"{idx}\n")

print(f"Generated ordered list of {len(ordered_frames)} frames")