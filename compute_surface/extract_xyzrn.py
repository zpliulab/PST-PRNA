from Bio.PDB import *
from default_config.dir_options import dir_opts
import os
radii = {}
radii["N"] = "1.540000"
radii["O"] = "1.400000"
radii["C"] = "1.740000"
radii["H"] = "1.200000"
radii["S"] = "1.800000"
radii["P"] = "1.800000"
radii["Z"] = "1.39"
radii["X"] = "0.770000"  ## Radii of CB or CA in disembodied case.
polarHydrogens = {}
polarHydrogens["ALA"] = ["H"]
polarHydrogens["GLY"] = ["H"]
polarHydrogens["SER"] = ["H", "HG"]
polarHydrogens["THR"] = ["H", "HG1"]
polarHydrogens["LEU"] = ["H"]
polarHydrogens["ILE"] = ["H"]
polarHydrogens["VAL"] = ["H"]
polarHydrogens["ASN"] = ["H", "HD21", "HD22"]
polarHydrogens["GLN"] = ["H", "HE21", "HE22"]
polarHydrogens["ARG"] = ["H", "HH11", "HH12", "HH21", "HH22", "HE"]
polarHydrogens["HIS"] = ["H", "HD1", "HE2"]
polarHydrogens["TRP"] = ["H", "HE1"]
polarHydrogens["PHE"] = ["H"]
polarHydrogens["TYR"] = ["H", "HH"]
polarHydrogens["GLU"] = ["H"]
polarHydrogens["ASP"] = ["H"]
polarHydrogens["LYS"] = ["H", "HZ1", "HZ2", "HZ3"]
polarHydrogens["PRO"] = []
polarHydrogens["CYS"] = ["H"]
polarHydrogens["MET"] = ["H"]



"""
xyzrn.py: Read a pdb file and output it is in xyzrn for use in MSMS
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

def extract_xyzrn(pdb_id, chain):
    """
        pdbfilename: input pdb filename
        xyzrnfilename: output in xyzrn format.
    """
    pdbfilename = os.path.join(dir_opts['protonated_pdb_dir'],pdb_id+'_'+chain+'.pdb')
    # pdbfilename = os.path.join(dir_opts['chain_pdb_dir'],pdb_id+'_'+chain+'.pdb')

    xyzrnfilename = os.path.join(dir_opts['xyzrn_dir'],pdb_id+'_'+ chain+'.xyzrn')
    if not os.path.exists(dir_opts['xyzrn_dir']):
        os.makedirs(dir_opts['xyzrn_dir'])

    parser = PDBParser()
    struct = parser.get_structure(pdbfilename, pdbfilename)

    out_list=[]
    for atom in struct.get_atoms():
        name = atom.get_name()
        residue = atom.get_parent()
        # Ignore hetatms.
        if residue.get_id()[0] != " ":
            continue
        resname = residue.get_resname().strip()
        chain = residue.get_parent().get_id()
        atomtype = name[0]

        coords=None
        if atomtype in radii:

            coords = "{:.06f} {:.06f} {:.06f}".format(
                atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
            )
            insertion="x"
            if residue.get_id()[2] != " ":
                insertion=residue.get_id()[2]
            full_id = "{}_{:d}_{}_{}_{}".format(
                chain, residue.get_id()[1], insertion, resname, name
            )
        if coords is not None:
            out_list.append(coords + " " + radii[atomtype] + " 1 " + full_id + "\n")
    outfile = open(xyzrnfilename, "w")
    outfile.writelines(out_list)
    outfile.close()
    #return xyzrnfilename