import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# --- Load the CP2K coordinates ---
elements = []
coords = []
with open('coord.inp','r') as f:
    inside = False
    for line in f:
        line = line.strip()
        if line.startswith('&COORD'):
            inside = True
            continue
        if line.startswith('&END COORD'):
            break
        if inside:
            tokens = line.split()
            elem = tokens[0]
            x, y, z = map(float, tokens[1:4])
            elements.append(elem)
            coords.append([x,y,z])
coords = np.array(coords)
elements = np.array(elements)

# --- Separate Uranium and Chlorine ---
mask_U  = elements == 'U'
mask_Cl = elements == 'Cl'
pts_U  = coords[mask_U]
pts_Cl = coords[mask_Cl]

# --- Plot in 3D ---
fig = plt.figure()
ax  = fig.add_subplot(111, projection='3d')

# Uranium atoms as circles
ax.scatter(pts_U[:,0], pts_U[:,1], pts_U[:,2],
           marker='o', label='U')

# Chlorine atoms as triangles
ax.scatter(pts_Cl[:,0], pts_Cl[:,1], pts_Cl[:,2],
           marker='^', label='Cl')

ax.set_xlabel('X (Å)')
ax.set_ylabel('Y (Å)')
ax.set_zlabel('Z (Å)')
ax.legend()
plt.tight_layout()
plt.show()
