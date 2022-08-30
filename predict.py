import os, re, random
import torch
import numpy as np
import pandas as pd
import torch.nn.functional as F

from nn_module.resnet18 import resnet18
from default_config.dir_options import dir_opts
from nn_module.function import f1, auc, save_ckpt, WeightedMultilabel,adjust_learning_rate
from sklearn.metrics import roc_curve, auc, recall_score, matthews_corrcoef,f1_score
import shutil
import torchvision.transforms as transform
import _pickle as cPickle

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print('Model is working on device', device)

class data_loader():
    def __init__(self, dataset, fold, is_train=True):
        self.is_train = is_train
        self.protein_list = []
        with open(dataset, 'r') as pid:
            for line in pid.readlines():
                self.protein_list.append(line.strip())

        self.data_all = []
        self.res_id_surface =[]
        self.res_id_all = []
        self.pdbid_all = []
        success = 0
        for pdb in self.protein_list:
            pdb=pdb.split('_')
            protein = pdb[0]
            chain = pdb[1] if len(pdb)==2 else ''
            pdb = protein +'_' + chain
            try:
                with open(os.path.join(dir_opts['data_label'], pdb), 'rb') as pid:
                    data0 = cPickle.load(pid)
                    data = data0['data'].transpose([0,3,1,2])
                    self.data_all.append(torch.tensor(data, dtype=torch.float32))
                    self.res_id_surface.append(data0['surface_res_id'])
                    self.res_id_all.append(data0['all_res_id'])
                    self.pdbid_all.append(protein)

            except:
                print(pdb,'not loaded')


    def __getitem__(self, item):
        return self.data_all[item], None, self.res_id_surface[item], self.res_id_all[item], self.pdbid_all[item]

    def __len__(self):
        return len(self.data_all)


def val_batch(protein_data, model, criterion):
    with torch.no_grad():
        model.eval()
        if torch.cuda.is_available():
            protein_data = protein_data.cuda()
        out = model(protein_data)

        return F.softmax(out, dim=1)

def predict():

    # model_dir = os.path.join(dir_opts['save_dir'], 'best_0.ckpt')
    model = resnet18()
    best_w = torch.load(dir_opts['model_dir'])
    model.load_state_dict(best_w['state_dict'])
    if torch.cuda.is_available():
        model = model.cuda()
    val_set = data_loader(dir_opts['PDB_list_to_predict'], fold=None,is_train=False)


    for feature, _, res_id_surface, res_id_all, pdb_id in val_set:
        criterion = WeightedMultilabel()
        proba_surface = val_batch(feature, model, criterion)
        proba_surface = proba_surface.detach().cpu().numpy()[:,1].tolist()
        if not os.path.exists(dir_opts['predict_dir']):
            os.mkdir(dir_opts['predict_dir'])

        proba_all = [0.0]*len(res_id_all)
        for index, res_id in enumerate(res_id_all):
            if res_id in res_id_surface:
                index_surface = res_id_surface.index(res_id)
                proba_all[index] = proba_surface[index_surface]
        data = {'res_id': res_id_all, 'RNA_binding_score': proba_all}
        df=pd.DataFrame(data)
        df.to_csv(os.path.join(dir_opts['predict_dir'], pdb_id+'.csv'))
        print(pdb_id,'predicted over')
        pass


if __name__ == '__main__':
    predict()
