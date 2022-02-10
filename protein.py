import os
import numpy as np
from subprocess import Popen, PIPE
from scipy.interpolate import griddata
import _pickle as cPickle
from sklearn.neighbors import KDTree

from Bio.PDB import *
from Bio.SeqUtils import seq1, seq3
from Bio.Seq import Seq
from Bio import SeqIO

from default_config.bin_path import bin_path
from default_config.dir_options import dir_opts

from compute_surface.protonate import protonate
from compute_surface.extractPDB import extractPDB
from compute_surface.extractRNA import extractcom
from compute_surface.extract_xyzrn import extract_xyzrn
from compute_surface.apply_msms import computeMSMS

from features.Hydrophobicity import kd_scale
from features.pKa import pKa_scale
from features.macc import macc
from features.Physicochemical import li_scale

class RBP():

    def __init__(self, protein_name):
        self.pdb_id, self.chain = protein_name.split('_') #complex, protein, rna

        #download pdb, call Reduce and MSMS
        self._download_pdb()

        self.model = self._load_pdb()

        self.RNA_space, self.RNA_chain = self._get_RNA_space()

        self.vertices, self.vert_info = self._get_surface() #ndarray (n*3), list(n) A_19_x_VAL_HG11

        #Extract the sequence of protien
        self.seq, self.index2resid, self.resid2Residue = self._get_pdb_seq(self.chain, kd_scale)

        #Get the coordinates of atom of all RNAs in complex
        self.geometric_center = self._get_geometric_center()

        #Get res_id on the surface
        self.surface_res = self._get_surface_res()

        #Get surface center of each residue to calculate the sampling density on the surface.
        self.res_to_vertice = self._get_res_to_vertice()

        #Calculate the label of each surface_res
        self.label, self.all_label = self._get_label()



    def _download_pdb(self):
        if not os.path.exists(dir_opts['raw_pdb_dir']):
            os.makedirs(dir_opts['raw_pdb_dir'])
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(self.pdb_id, pdir=dir_opts['raw_pdb_dir'], file_format='mmCif')

    def _get_surface(self):
        extractPDB(self.pdb_id + '_' + self.chain)
        protonate(self.pdb_id+'_'+self.chain)
        extract_xyzrn(self.pdb_id, self.chain)
        vertices, faces1, normalv1, vert_info = computeMSMS(self.pdb_id, self.chain)
        return vertices, vert_info

    def _load_pdb(self):
        pdb_file = os.path.join(dir_opts['raw_pdb_dir'], self.pdb_id+'.cif')
        parser = MMCIFParser(QUIET=True)
        struct = parser.get_structure(self.pdb_id, pdb_file)
        model = Selection.unfold_entities(struct, "M")[0]
        return model
    
    def _get_pdb_seq(self, chain, scale):
        chain = self.model.child_dict[chain]
        res_seq = ''
        index2resid = {}
        resid2structure = {}
        index = 0
        for Residue in chain.child_list:
            if Residue.get_resname() in scale.keys():
                res_seq += (seq1(Residue.get_resname()))
                index2resid[index] = Residue.get_id()[1]
                resid2structure[Residue.get_id()[1]] = Residue
                index += 1
        return res_seq, index2resid, resid2structure

    def _get_RNA_space(self):
        RNA_list = ["A", "C", "G", "U"]
        RNA_chain = set()
        atom_list = []
        for chain in self.model.child_list:
            if chain.id == self.chain:
                continue
            for res in chain.child_list:
                atom_type = res.resname.strip()
                if atom_type in RNA_list:
                    for atom in res:
                        atom_list.append(atom.coord)
                        RNA_chain.add(chain.id)
        chain = '_'.join(RNA_chain)
        return atom_list, chain

    def _get_label(self):
        rna_tree = KDTree(self.RNA_space)
        label_dict = {}
        for res_id in self.surface_res:
            res = self.resid2Residue[res_id]
            res_coord = []
            for atom in res.child_list:
                res_coord.append(atom.coord)
            d, t = rna_tree.query(res_coord)
            if np.min(d) < 5.0:
                label_dict[res_id] = 1
            else:
                label_dict[res_id] = 0
        label_all = []
        for index in range(len(self.seq)):
            if self.index2resid[index] in self.surface_res:
                label_all.append(label_dict[self.index2resid[index]])
            else:
                label_all.append(0)
        return label_dict, np.array(label_all)

    def _get_surface_res(self):
        surface_res = set()
        for item in self.vert_info:
            resid = int(item.split('_')[1])
            surface_res.add(resid)
        return list(surface_res)

    def _get_res_to_vertice(self):
        res_to_vertices = {}
        for index, item in enumerate(self.vert_info):
            res_id = int(item.split('_')[1])
            if res_id in res_to_vertices.keys():
                res_to_vertices[res_id].append(self.vertices[index,:])
            else:
                res_to_vertices[res_id] = [self.vertices[index,:],]
        res_to_vertice = {}
        for item in res_to_vertices.keys():
            res_surface_array = np.array(res_to_vertices[item])
            res_to_vertice[item] = np.mean(res_surface_array,axis=0)
        return res_to_vertice

    def _get_geometric_center(self):
        return self.vertices.mean(axis=0)

    def _call_dssp(self):
        dssp_bin = bin_path['DSSP']
        pdb_file = os.path.join(dir_opts['chain_pdb_dir'], self.pdb_id + '_'+self.chain + '.pdb')
        if not os.path.exists(dir_opts['dssp']):
            os.mkdir(dir_opts['dssp'])
        out_file = os.path.join(dir_opts['dssp'], self.pdb_id+'_'+self.chain+'.dssp')

        args = [dssp_bin, '-i', pdb_file, '-o', out_file]
        p2 = Popen(args, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p2.communicate()
        rASA = {}
        secondary = {}
        flag_1 = False
        flag_2 = False
        with open(out_file) as pid:
            for line in pid.readlines():
                if line[2] == '#':
                    flag_1 = True
                    continue
                if flag_1 and line[11] == self.chain[-1]:
                    flag_2 = True
                    res_name = line[13]
                    C_list = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
                    if res_name in C_list:
                        res_name = 'C'
                    if res_name == 'X':
                        continue
                    res_id = int(line[5:10])
                    rASA[res_id] = float(line[34:39]) / float(macc[res_name])
                    if line[16] in ('G', 'H', 'I'):
                        secondary[res_id] = [1, 0, 0]
                    if line[16] in ('E', 'B'):
                        secondary[res_id] = [0, 1, 0]
                    if line[16] in ('T', 'S', ' '):
                        secondary[res_id] = [0, 0, 1]
                if flag_1 and flag_2 and line[13:15] == '!*':
                    break
        return rASA, secondary

    def _call_psiblast(self):
        psiblast_bin = bin_path['PSIBLAST']
        uniprot_database = bin_path['PSIBLAST_DATABASE']

        if not os.path.exists(dir_opts['blast_fasta_dir']):
            os.mkdir(dir_opts['blast_fasta_dir'])
        if not os.path.exists(dir_opts['blast_pssm_dir']):
            os.mkdir(dir_opts['blast_pssm_dir'])

        #save fasta file
        seq = Seq(self.seq)
        seq_io = SeqIO.SeqRecord(seq, name = self.pdb_id + '_' + self.chain)
        fasta_dir = os.path.join(dir_opts['blast_fasta_dir'], self.pdb_id + '_' + self.chain + '.fasta')
        if not os.path.exists(fasta_dir):
            SeqIO.write(seq_io, fasta_dir, "fasta")

        #make pssm_file
        pssm_dir = os.path.join(dir_opts['blast_pssm_dir'],self.pdb_id+'_'+self.chain+'.pssm')
        args = [psiblast_bin, "-db", uniprot_database, "-query", fasta_dir, "-evalue", '0.001', "-num_iterations", '3', "-out_ascii_pssm", pssm_dir]
        if not os.path.exists(pssm_dir):
            p2 = Popen(args, stdout=PIPE, stderr=PIPE)
            print('calling psiblast start')
            stdout, stderr = p2.communicate()
            print('calling psiblast over')

        #sparse pssm_file
        RNA_order = 'ARNDCQEGHILKMFPSTWYV'
        pssm_feature = {}
        with open(pssm_dir, 'r') as pid:
            for index, line in enumerate(pid.readlines()[3:]):
                if line == '\n':
                    break
                pssm_acid, feature = self._sparse_pssm_line(line)
                # amino_acid_position = RNA_order.index(pssm_acid)
                pssm_feature[self.index2resid[index]] = feature
        return pssm_feature

    def _call_hhm(self):
        HHblits = bin_path['HHblits']
        HHblits_DB = bin_path['HHblits_DB']

        if not os.path.exists(dir_opts['blast_fasta_dir']):
            os.mkdir(dir_opts['blast_fasta_dir'])
        if not os.path.exists(dir_opts['hhm_dir']):
            os.mkdir(dir_opts['hhm_dir'])

        #save fasta file
        seq = Seq(self.seq)
        seq_io = SeqIO.SeqRecord(seq, name = self.pdb_id + '_' + self.chain)
        fasta_dir = os.path.join(dir_opts['blast_fasta_dir'], self.pdb_id + '_' + self.chain + '.fasta')
        if not os.path.exists(fasta_dir):
            SeqIO.write(seq_io, fasta_dir, "fasta")

        hhm_dir = os.path.join(dir_opts['hhm_dir'], self.pdb_id+'_'+self.chain+'.hhm')

        args = [HHblits, '-d', HHblits_DB, '-i', fasta_dir, '-ohhm', hhm_dir]
        if not os.path.exists(hhm_dir):
            p2 = Popen(args, stdout=PIPE, stderr=PIPE)
            print('calling hhblits start')
            stdout, stderr = p2.communicate()
            print('calling hhblits over')
        hhm_feature = {}
        hhm_dir = os.path.join(dir_opts['hhm_dir'], self.pdb_id+'_'+self.chain+'.hhm')

        with open(hhm_dir, 'r') as f:
            text=f.readlines()
        hhm_begin_line=0
        hhm_end_line=0
        for i in range(len(text)):
            if '#' in text[i]:
                hhm_begin_line=i + 5
            elif '//' in text[i]:
                hhm_end_line=i
        hhm=np.zeros([int((hhm_end_line - hhm_begin_line) / 3), 30])

        axis_x=0
        for i in range(hhm_begin_line, hhm_end_line, 3):
            line1=text[i].split()[2:-1]
            line2=text[i + 1].split()
            axis_y=0
            for j in line1:
                if j == '*':
                    hhm[axis_x][axis_y]=9999 / 10000.0
                else:
                    hhm[axis_x][axis_y]=float(j) / 10000.0
                axis_y+=1
            for j in line2:
                if j == '*':
                    hhm[axis_x][axis_y]=9999 / 10000.0
                else:
                    hhm[axis_x][axis_y]=float(j) / 10000.0
                axis_y+=1
            axis_x+=1
        hhm=(hhm - np.min(hhm)) / (np.max(hhm) - np.min(hhm))
        for index in range(int((hhm_end_line - hhm_begin_line) / 3)):
            resid = self.index2resid[index]
            hhm_feature[resid] = hhm[index,:].tolist()

        return hhm_feature

    def _sparse_pssm_line(self, line):
        a = line.split()
        pssm_str = a[2:22]
        pssm_value = []
        for item in pssm_str:
            # pssm_value.append(1/(1+np.exp(float(item))))
            pssm_value.append(float(item))
        amino_acid = a[1]
        return amino_acid, pssm_value

    def _get_pKa_feature(self):
        #used to calculate pKa feature
        pKa_mapping = {}
        for item in pKa_scale.keys():
            pKa_mapping[item] = ((pKa_scale[item][0] - 1.82) / (2.38-1.82), (pKa_scale[item][1] - 8.80) / (10.96 - 8.80),
                                 (pKa_scale[item][2] - 3.65) / (12.48-3.65))
        return pKa_mapping

    def _get_Physicochemical_feature(self):
        physio_scale = {}
        for item in li_scale.keys():
            physio_scale[seq3(item).upper()] = ((li_scale[item][0] - (-1))/(1-(-1)), ((li_scale[item][1] - 2)/(4-2)))
        return physio_scale

    def _get_hydrophobicity_feature(self):
        #used to calculate hydrophily
        hydrophobicity_scale = {}
        for item in kd_scale.keys():
            hydrophobicity_scale[item] = (kd_scale[item] - (-4.5)) /(4.5 - (-4.5))
        return hydrophobicity_scale

    def _get_fea_mapping(self):
        #call dssp and psiblast
        self.feature_rASA, self.feature_secondary = self._call_dssp()
        self.feature_pssm = self._call_psiblast()
        self.feature_hhm = self._call_hhm()
        self.feature_pKa = self._get_pKa_feature() #shape: (length_seq, 3)
        self.feature_physicochemical = self._get_Physicochemical_feature() #(length_seq, 2)
        self.feature_hydophobicity = self._get_hydrophobicity_feature() #(length_seq, 1)

    def _get_rotation_matrix(self, res_id):
        vertex_center = self.res_to_vertice[res_id]
        vertex_center = vertex_center - self.geometric_center
        rotation_x = np.arctan2(vertex_center[2], vertex_center[1])
        if rotation_x >= 0:
            rotation_x = np.pi * 2 - rotation_x
        else:
            rotation_x = - rotation_x
        rotation_matrix_x = np.array([[1, 0, 0],
                                    [0, np.cos(rotation_x), np.sin(rotation_x)],
                                    [0, -np.sin(rotation_x), np.cos(rotation_x)],
                                    ])
        vertex_center_rotation_x = np.matmul(vertex_center, rotation_matrix_x)
        rotation_z = np.arctan2(vertex_center_rotation_x[1], vertex_center_rotation_x[0])
        if rotation_z >= 0:
            rotation_z = np.pi * 2 - rotation_z
        else:
            rotation_z = -rotation_z
        rotation_matrix_z = np.array([[np.cos(rotation_z), np.sin(rotation_z), 0],
                                    [-np.sin(rotation_z), np.cos(rotation_z), 0],
                                    [0, 0, 1],
                                    ])
        return rotation_matrix_x, rotation_matrix_z

    def _get_polar_coord(self, res_id):
        rotation_x, rotation_z = self._get_rotation_matrix(res_id)
        coord = self.vertices - self.geometric_center
        coord_rotation_x = np.matmul(coord,rotation_x)
        coord = np.matmul(coord_rotation_x, rotation_z)
        #get all lat lon
        x, y, z, radius = coord[:,0], coord[:,1], coord[:,2], np.sqrt(np.sum((coord ** 2), axis=1))
        lat = np.arcsin(z / radius) * 180 / np.pi  # Latitudes
        lon = np.arctan2(y, x) * 180 / np.pi  # Longitudes
        lon_lat = np.hstack([lon[:, None], lat[:, None]])
        #get min and max lat and lon
        center_point = self.res_to_vertice[res_id] - self.geometric_center
        center_point = np.matmul(center_point, rotation_x)
        center_point = np.matmul(center_point, rotation_z)
        kdtree = KDTree(coord)
        (ind, dist) = kdtree.query_radius(center_point[None,:], 12, return_distance = True)
        coord = coord[ind.tolist()[0],:]
        x, y, z, radius = coord[:,0], coord[:,1], coord[:,2], np.sqrt(np.sum((coord ** 2), axis=1))
        lat = np.arcsin(z / radius) * 180 / np.pi  # Latitudes
        lon = np.arctan2(y, x) * 180 / np.pi  # Longitudes
        min_max_lon_lat = (np.min(lon), np.max(lon), np.min(lat), np.max(lat))
        return lon_lat, min_max_lon_lat

    def _interpolation(self, polar_pix, min_max_lon_lat, feature):
        x_lons = np.linspace(min_max_lon_lat[0], min_max_lon_lat[1], 32)
        y_lats = np.linspace(min_max_lon_lat[2], min_max_lon_lat[3], 32)
        MhpS = griddata(polar_pix, feature, (x_lons[None, :], y_lats[:, None]), method='nearest')
        return MhpS

    def _get_graph(self):   
        self._get_fea_mapping()
        coord = []
        graph = []
        for index, vert in enumerate(self.vert_info):
            try:
                res_name = vert.split('_')[3]
                res_id = int(vert.split('_')[1])
                fea = []
                fea.append(self.feature_rASA[res_id]) #1
                fea.extend(self.feature_secondary[res_id]) #3
                fea.extend(self.feature_pssm[res_id])#20
                fea.extend(self.feature_hhm[res_id])#30
                fea.extend(self.feature_pKa[res_name])#3
                fea.extend(self.feature_physicochemical[res_name])#2
                fea.append(self.feature_hydophobicity[res_name])#1
                coord.append(self.vertices[index, :])
                graph.append(fea)
            except:
                pass
        self.vertices = coord
        return np.array(graph, dtype=np.float)

    def get_data(self):
        self.graph = self._get_graph()
        data = []
        label = []
        res_surface = []
        for index, res_id in enumerate(self.surface_res):
            try:
                lon_lat, min_max_lon_lat = self._get_polar_coord(res_id)
                features = self._interpolation(lon_lat, min_max_lon_lat, self.graph)
                data.append(features)
                label.append(self.label[res_id])
                res_surface.append(res_id)
            except:
                print(self.pdb_id + '_' + self.chain, index, res_id, 'wrong')
                continue
        batch_dict = {'data':np.array(data),  'label':np.array(label), 'res_id':res_surface}
        with open(os.path.join(dir_opts['data_label'], self.pdb_id + '_' + self.chain ), 'wb') as pid:
            cPickle.dump(batch_dict, pid)