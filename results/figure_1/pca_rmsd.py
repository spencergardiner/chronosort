import numpy as np
import matplotlib.pyplot as plt

# Load the RMSD data from the .xvg file
data = np.loadtxt('rmsd.xvg', comments=['@', '#'])

# Extract time and RMSD columns (assuming 2 columns: time, RMSD)
time = data[:, 0]  # First column is time (ns)
rmsd = data[:, 1]  # Second column is RMSD (nm)

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(time, rmsd, label='RMSD', color='limegreen')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (nm)')
plt.title('RMSD vs Time')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save the figure
plt.savefig('rmsd_plot.png', dpi=600)
plt.show()
