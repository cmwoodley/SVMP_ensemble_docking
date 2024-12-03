from pymol import cmd
from pymol.cmd import util
import os
from openbabel import pybel

def replace_v_with_zn(pdb_filepath, output_filepath):
    with open(pdb_filepath, 'r') as file:
        lines = file.readlines()

    with open(output_filepath, 'w') as file:
        for line in lines:
            # Check if the line corresponds to the atom type "V"
            if line.startswith("HETATM") and "V" in line[76:78]:
                # Replace the atom name and element symbol from " V " to "ZN "
                line = line[:12] + " ZN " + line[16:76] + "ZN" + line[78:]
            file.write(line)

def main():
    for prot_name in ["PI", "PII", "PIII_1", "PIII_2", "PIII_3"]:
        prot_dir = f"../api_docking/{prot_name}/"
        for lig_name in ["DC174", "XL784","prinomastat"]:
            best_dir = f"../api_docking/{prot_name}/{lig_name}/best_poses/"

            for i,file in enumerate([x for x in os.listdir(best_dir) if x.endswith(".pdb")]):
                filename = os.path.join(best_dir, file)

                # Read the molecule from the PDB file
                mol = next(pybel.readfile("pdb", filename))

                # Filter out atoms with type "Du"
                atoms_to_remove = [atom for atom in mol.atoms if atom.type == "Du"]

                # Remove the "Du" atoms
                for atom in atoms_to_remove:
                    mol.OBMol.DeleteAtom(atom.OBAtom)

                mol.write("pdb", filename, overwrite=True)  

                replace_v_with_zn(filename, filename) 
                
                cmd.load(filename, f"{lig_name}_pose_{i}")

            cmd.group(lig_name, " ".join([f"{lig_name}_pose_{i}" for i in range(10)]))

        cmd.select("want","polymer or name ZN or name CA or organic")
        cmd.select("not_want", "not want")
        cmd.select("close_res","byres polymer near_to 3.0 of organic")

        cmd.remove("not_want")

        cmd.bg_color("white")
        cmd.set("cartoon_color", "gray90")
        cmd.show("sticks","close_res")
        cmd.color_deep("gray90", 'close_res', 0)
        util.cnc("close_res",_self=cmd)
        cmd.hide("sticks","hydrogen and (elem C extend 1)")
        cmd.delete("not_want")
        cmd.delete("want")
        cmd.delete("close_res")

        cmd.save(os.path.join(prot_dir, f"{prot_name}.pse"))
        cmd.delete("all")

if __name__ == "__main__":
    main()