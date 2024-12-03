from ccdc.docking import Docker
from ccdc.io import MoleculeReader, MoleculeWriter
import os

# Define directories and constraint atoms
dock_dir = "../api_docking/"
constraint_atoms = {
    "DC174": [43, 46, 47],
    "XL784": [14, 17, 18],
    "prinomastat": [20,23,24]
}

bs_atoms = {
    "PI":2208,
    "PII":2496,
    "PIII_1":2590,
    "PIII_2":1162,
    "PIII_3":2244,
}

protein_glu = {
    "PI":2223,
    "PII":2182,
    "PIII_1":2243,
    "PIII_2":826,
    "PIII_3":2260,    
}

def main():
    contents = os.listdir()

    for prot_name in ["PI", "PII", "PIII_1", "PIII_2", "PIII_3"]:
        for lig_name in ["DC174", "XL784", "prinomastat"]:

            print(f"Docking {lig_name} into {prot_name}")

            docker = Docker()
            settings = docker.settings

            ligand_file = os.path.join(dock_dir, "ligands", f"{lig_name}.sdf")
            prot_dir = os.path.join("../data", prot_name, "clust_outputs")
            outdir = os.path.join(dock_dir, prot_name, lig_name)
            os.makedirs(outdir, exist_ok=True)

            ligand = MoleculeReader(ligand_file)[0]
            settings.add_ligand_file(ligand_file, 20)

            protein_files = [x for x in os.listdir(prot_dir) if x.endswith(".pdb")]
            print(prot_dir)
            for i, protein_file in enumerate(protein_files): 
                settings.add_protein_file(os.path.join("../data", prot_name, "clust_outputs", protein_file))
                protein = settings.proteins[i]

                zn = next(j for j,a in enumerate(protein.atoms) if a.residue_label == "ZN2")
                c0 = settings.DistanceConstraint(settings.proteins[i].atoms[zn], ligand.atoms[constraint_atoms[lig_name][0]], (1.5, 2.5), 10)            
                c1 = settings.DistanceConstraint(settings.proteins[i].atoms[zn], ligand.atoms[constraint_atoms[lig_name][1]], (1.5, 2.5), 10)     
                c2 = settings.DistanceConstraint(settings.proteins[i].atoms[protein_glu[prot_name]], ligand.atoms[constraint_atoms[lig_name][2]], (1.5, 2.5), 10)     
                settings.protein_files[i].add_constraint(c0)       
                settings.protein_files[i].add_constraint(c1)       
                settings.protein_files[i].add_constraint(c2)       

            settings.binding_site = settings.BindingSiteFromAtom(settings.proteins[0], settings.proteins[0].atoms[bs_atoms[prot_name]], 6)
            settings.fitness_function = 'plp'
            settings.autoscale = 10.0
            settings.early_termination = False
            settings.write_options = 'NO_LOG_FILES'
            settings.output_directory = outdir
            # settings.output_file = 'poses.sdf'
            settings.flip_amide_bonds = True
            settings.flip_free_corners = True
            settings.flip_pyramidal_nitrogen = True
            settings.save_lone_pairs = True

            results = docker.dock()
            docked_ligands = results.ligands
            pdb_dir = os.path.join(outdir, "best_poses")
            os.makedirs(pdb_dir, exist_ok=True)
            for i in range(10):
                complex = results.make_complex(docked_ligands[i])
                with MoleculeWriter(os.path.join(pdb_dir, f"pose_{i+1}.pdb")) as mol_writer:
                    mol_writer.write(complex)

            # Handle api_gold.conf file
            conf_path = os.path.join(outdir, "api_gold.conf")
            if os.path.exists(conf_path):
                os.remove(conf_path)
            os.rename("api_gold.conf", conf_path)

            # Cleanup any new files created during the process
            for file in os.listdir():
                if file not in contents:
                    os.remove(file)

if __name__ == "__main__":
    main()
