import os
from default_config.dir_options import dir_opts
from Bio.PDB import *
from features.Hydrophobicity import kd_scale
RNA_list=["A", "C", "G", "U"]

#RNA_chains are only used to transform from cif format to pdb format
RNA_chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"  or atom.get_altloc() == "1"
# Exclude disordered atoms.
def extractPDB(pair):

    pairs = pair.split('_')
    pdb_id = pairs[0]
    protein_chain = pairs[1]
    infilename = os.path.join(dir_opts['raw_pdb_dir'],pdb_id+'.cif')
    if not os.path.exists(dir_opts['chain_pdb_dir']):
        os.mkdir(dir_opts['chain_pdb_dir'])

    outfilename = os.path.join(dir_opts['chain_pdb_dir'], pair+'.pdb')
    # extract the chain_ids from infilename and save in outfilename.

    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure(infilename, infilename)
    model = Selection.unfold_entities(struct, "M")[0]
    # Select residues to extract and build new structure
    structBuild = StructureBuilder.StructureBuilder()
    structBuild.init_structure("output")
    structBuild.init_seg(" ")
    structBuild.init_model(0)
    outputStruct = structBuild.get_structure()
    RNA_index = 0
    for chain in model:
        if chain.get_id()==protein_chain:
            structBuild.init_chain(chain.get_id()[-1])
            for residue in chain:
                if residue.get_resname().upper() in kd_scale.keys():
                    outputStruct[0][chain.get_id()[-1]].add(residue)

    # Output the selected residues
    pdbio = PDBIO()
    pdbio.set_structure(outputStruct)
    pid = open(outfilename,'w')
    pdbio.save(pid, select=NotDisordered())
    pid.close()