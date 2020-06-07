from pathlib import Path
from typing import List

def pymol_script_writer(out_file: str, pdb_id: str, pos: list): 
    path = Path(out_file)

    with open(path, "w") as file: 
        file.write(f"fetch {pdb_id}\n\n")

        for i,_ in enumerate(pos): 
            file.write(f"create mot{pos[i]}, resi {pos[i]}-{pos[i]+4} \n")
            
        file.write("\nhide_all\n")
        file.write("\n")

        for i,_ in enumerate(pos): 
            file.write(f"show cartoon, resi{pos[i]}\n")
            
        file.write("\n")
        for i,_ in enumerate(pos): 
            file.write(f"color red, mot{pos[i]}\n")