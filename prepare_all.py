import os, sys
from protein import RBP
import numpy as np
from default_config.dir_options import dir_opts
if __name__ == "__main__":
    dataset_path = 'data/pdbid_chain/RBP03.txt'
    #dataset_path = 'data/pdbid_chain/RBP407_09.txt'

    protein_list = []
    with open(dataset_path, 'r') as pid:
        for line in pid.readlines():
            protein_list.append(line.strip())

    if not os.path.exists(dir_opts['data_label']):
        os.mkdir(dir_opts['data_label'])

    # pdb_slurm = sys.argv[1].strip()
    print(len(protein_list))
    for index, item in enumerate(protein_list):
        # if item == pdb_slurm:
            if os.path.exists(os.path.join(dir_opts['data_label'], item)):
                continue
            print(index, item)
            rbp = RBP(item)
            rbp.get_data()