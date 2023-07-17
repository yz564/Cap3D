import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def read_csv_with_multiple_tables(file_path):
    table_dfs = []
    table_names = []

    with open(file_path, 'r') as file:
        table_started = False
        table_name = ""
        table_content = ""

        for line in file:
            line = line.strip()

            if line:
                if not table_started:
                    table_name = line
                    table_started = True
                else:
                    table_content += line + "\n"
            elif table_started:
                table_dfs.append(pd.read_csv(io.StringIO(table_content.strip())))
                table_names.append(table_name)
                table_started = False
                table_content = ""
                table_name=""
        if (table_name): #add the last table if the last line is not \n
                table_dfs.append(pd.read_csv(io.StringIO(table_content.strip())))
                table_names.append(table_name)
        
    return table_dfs, table_names
    
def find_table_by_name(tables, table_names, target_table_name):
    for i, name in enumerate(table_names):
        if name == target_table_name:
            return tables[i]
    return None
    
def plot_triangle_mesh(vertices, triangles, values):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot vertices
    ax.scatter3D(vertices[:, 0], vertices[:, 1], vertices[:, 2], s=1, c='black')

    # Plot triangles
    mesh = [vertices[tri] for tri in triangles]
    tri_collection = Poly3DCollection(mesh, linewidths=0.5, edgecolors='black', alpha=1, cmap='seismic')
    tri_collection.set_array(values)
    ax.add_collection3d(tri_collection)
    
    # Set axis ratio to be the same (equal aspect ratio)
    #ax.set_box_aspect((np.ptp(vertices[:, 0]), np.ptp(vertices[:, 1]), np.ptp(vertices[:, 2])))

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Charge density')
    
    cbar = plt.colorbar(tri_collection)
    cbar.set_label('C/m^2', rotation=270, labelpad=15)

    plt.savefig('charge_density.png')
    plt.show()

if __name__ == "__main__":
    csv_file = "mesh_charge_density.csv"  # Replace this with your CSV file path
    tables, table_names = read_csv_with_multiple_tables(csv_file)
    vertices = find_table_by_name(tables, table_names, 'coordinates')
    triangles = find_table_by_name(tables, table_names, 'elements')
    charge = find_table_by_name(tables, table_names, 'charge_density')
    plot_triangle_mesh(vertices.values, triangles.values, charge.values.flatten())

