import mdtraj as md
import scipy
from sklearn.cluster import AgglomerativeClustering
import numpy as np 
import matplotlib.pyplot as plt
import os
import pymol
from pymol import cmd


def load_and_align_pdbs(folder_path):
    
    # Get all PDB files in the folder
    pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    
    if not pdb_files:
        return
    
    # Load the first PDB file to select atoms within 15 angstroms of ZN
    first_pdb = os.path.join(folder_path, pdb_files[0])
    print(first_pdb)
    cmd.load(first_pdb, "first")

    # Select atoms within 15 angstroms of ZN
    zn_selection = 'byres polymer near_to 12.0 of name ZN'
    cmd.select('selection', zn_selection)
    for pdb_file in pdb_files:
        pdb_path = os.path.join(folder_path, pdb_file)
        cmd.load(pdb_path, pdb_file)
        # Align current structure to the first structure based on the selection
        cmd.align(pdb_file, 'selection')
        
        # Save the aligned structure
        aligned_pdb_path = os.path.join(folder_path, pdb_file)
        cmd.save(aligned_pdb_path, pdb_file)
        
        # Delete the loaded structure to free memory
        cmd.delete(pdb_file)

    cmd.delete('all')

def main():
    for protein in os.listdir("../data/"):
        traj_dir = "../data/{}/raw_traj/".format(protein)
        if len(os.listdir(traj_dir)) != 0:
            for i, file in enumerate([x for x in os.listdir(traj_dir) if x.endswith("xtc")]):
                if i == 0:
                    t = md.load(traj_dir+file,top=traj_dir+"newbox.gro")
                    t = t[100:]
                if i != 0:
                    _t = md.load(traj_dir+file,top=traj_dir+"newbox.gro")
                    t= md.join([t,_t[100:]])

                    
            topology = t.topology

            # Select atoms within 15 angstroms of the zinc atom
            zn_atom = [atom.index for atom in topology.atoms if str(atom.residue) == "ZN2"][0]

            # Compute neighbors within 15 angstroms of the zinc atom in the first frame
            query_indices = [zn_atom]
            haystack_indices = topology.select('all')
            neighbors = md.compute_neighbors(t, cutoff=1.2, query_indices=query_indices, haystack_indices=haystack_indices, periodic=False)[0]

            # Create a new trajectory with only the selected atoms
            t_selected = t.atom_slice(neighbors)

            # Superpose the selected atoms trajectory on the first frame
            t_selected.superpose(t_selected, 0)
            # Calculate the RMSD distance matrix
            distances = np.empty((t_selected.n_frames, t_selected.n_frames))
            for i in range(t_selected.n_frames):
                distances[i] = md.rmsd(t_selected, t_selected, i)

            print('Max pairwise RMSD: %f nm' % np.max(distances))
            #Save cluster summary and bar chart

            clust = AgglomerativeClustering(distance_threshold=None, n_clusters=12, metric="precomputed", linkage="average")
            assignments = clust.fit_predict(distances)


            _, ax = plt.subplots(1,1,figsize=(6,4))
            scatter = ax.scatter(np.arange(len(distances[0])),distances[0], c=assignments)
            unique_assignments, _ = np.unique(assignments, return_inverse=True)
            legend_colors = scatter.to_rgba(np.unique(assignments))

            # Plot a single point for each unique assignment color with the corresponding label
            for i, (assignment, color) in enumerate(zip(unique_assignments, legend_colors)):
                scatter = ax.scatter([], [], c=[color], label=assignment)

            # Show the legend
            ax.set_xlabel("Concatentated Frame")
            ax.set_ylabel("RMSD (nm)")
            
            ax.legend(loc='right')
            plt.savefig("../reports/{}_rmsd.png".format(protein))
            centroids = []
            for j in range(len(np.unique(assignments))):
                distances_0 = distances[:,np.where(assignments == j)[0]]
                distances_0 = distances_0[np.where(assignments == j)[0]]
                beta = 1
                index = np.exp(-beta*distances_0 / distances_0.std()).sum(axis=1).argmax()
                centroids.append(index)

            centroids = [centroids[i]+min(np.where(assignments == i)[0]) for i in range(len(centroids))]

            for i in range(len(centroids)):
                filename = "../data/{}/clust_outputs/cluster_{}.pdb".format(protein,i)
                t[centroids[i]].save_pdb(filename)

            with open(filename,"r") as f:
                pdb = f.readlines()

            with open(filename,"w") as f:
                for i,line in enumerate(pdb):
                    if line.find("ATOM") != -1 and line.find("ZN") != -1 or line.find("HETATM") != -1 and line.find("ZN") != -1:
                        f.writelines(line[:76]+"ZN"+line[78:])
                    else:
                        f.writelines(line)

            load_and_align_pdbs(f"../data/{protein}/clust_outputs/")
            filename = "../reports/cluster_populations_{}.txt".format(protein)

            all_frames = np.stack(np.unique(assignments, return_counts=True)).T[:,1].sum()

            with open(filename,"w") as f:
                f.writelines("Summary of {} clustering\n".format(protein))
                for i,j in np.stack(np.unique(assignments, return_counts=True)).T:
                    f.writelines("Cluster_{}: {} frames ({:.1%})\n".format(i,j, j/all_frames))
        else:
            continue

if __name__ == "__main__":
    main()