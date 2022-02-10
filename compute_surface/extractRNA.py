import os
from default_config.dir_options import dir_opts
from Bio.PDB import *
from features.Hydrophobicity import kd_scale
from features.macc import macc
RNA_list=["A", "C", "G", "U"]

#RNA_chains are only used to transform from cif format to pdb format
RNA_chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"  or atom.get_altloc() == "1"
# Exclude disordered atoms.
def extractcom(pairs):
    pdb_id=pairs.split('_')[0]
    # extract the chain_ids from infilename and save in outfilename.
    infilename=os.path.join(dir_opts['raw_pdb_dir'], pdb_id + '.cif')
    parser=MMCIFParser(QUIET=True)
    struct=parser.get_structure(infilename, infilename)
    model = Selection.unfold_entities(struct, "M")[0]

    protein_chain=pairs.split('_')[1]


    for RNA in pairs.split('_')[2:]:
        pair = pdb_id + '_' + protein_chain +'_' + RNA
        if not os.path.exists(dir_opts['chain_pdb_dir']):
            os.mkdir(dir_opts['chain_pdb_dir'])

        outfilename = os.path.join(dir_opts['chain_pdb_dir'], pair +'.pdb')

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

            elif chain.get_id()==RNA:
                RNA_chain = RNA_chains[RNA_index]
                if RNA_chain == protein_chain[-1]:
                    RNA_index += 1
                    RNA_chain = RNA_chains[RNA_index]
                structBuild.init_chain(RNA_chain)
                for residue in chain:
                    if residue.get_resname() in RNA_list:
                        # residue.parent ={}
                        outputStruct[0][RNA_chain].add(residue)

                        # residue.id = (residue.id[0],RNA_res_id,residue.id[2])
                RNA_index +=1
        # Output the selected residues
        pdbio = PDBIO()
        pdbio.set_structure(outputStruct)
        pid = open(outfilename,'w')
        pdbio.save(pid, select=NotDisordered())
        pid.close()