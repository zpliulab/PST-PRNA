from default_config.dir_options import dir_opts
from default_config.bin_path import bin_path
from compute_surface.read_msms import read_msms
from subprocess import Popen, PIPE
import os
import random
def computeMSMS(pdb, chain):
    xyzrn = os.path.join(dir_opts['xyzrn_dir'], pdb+'_'+chain+'.xyzrn')
    msms_file_base = os.path.join(dir_opts['msms_dir'],  pdb+'_'+chain + str(random.randint(1,10000000)))
    if not os.path.exists(dir_opts['msms_dir']):
        os.makedirs(dir_opts['msms_dir'])

    msms_bin = bin_path['MSMS']
    # Now run MSMS on xyzrn file
    FNULL = open(os.devnull, 'w')
    args = [msms_bin, "-density", "3.0", "-hdensity", "3.0", "-probe",\
                    "1.5", "-if",xyzrn,"-of",msms_file_base, "-af", msms_file_base]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()
    vertices, faces, normalv, res_id = read_msms(msms_file_base)
    return vertices, faces, normalv, res_id

def get_asa(file_base, protein_chain:str):
    if '_' in protein_chain:
        protein_chain = protein_chain.split('_')[0]
    with open(file_base + '.area', 'r') as pid:
        lines = pid.readlines()
    sas = {}
    for line in lines[1:-1]:
        area = float(line[15:23])
        chain = line.strip().split()[-1].split('_')[0]
        res_id = int(line.strip().split()[-1].split('_')[1])
        if chain != protein_chain:
            continue
        if res_id not in sas.keys():
            sas[res_id] = area
        else:
            sas[res_id] += area
    return sas