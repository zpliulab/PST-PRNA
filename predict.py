import os, re
import torch
import numpy as np

import torch.nn.functional as F

from nn_module.resnet18 import resnet18
from default_config.dir_options import dir_opts
from nn_module.function import f1, auc, save_ckpt, WeightedMultilabel,adjust_learning_rate
from data_loader.data_loader import data_loader
from sklearn.metrics import roc_curve, auc, recall_score, matthews_corrcoef,f1_score
import shutil

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print('Model is working on device', device)

def val_batch(protein_data, site_label, model, criterion):
    with torch.no_grad():
        model.eval()
        if torch.cuda.is_available():
            protein_data = protein_data.cuda()
            site_label = site_label.cuda()
        out = model(protein_data)
        loss = criterion(out, site_label)
        loss = loss.data.item()
        return loss, F.softmax(out, dim=1)

def predict(dataset, output_base):

    model_dir = os.path.join(dir_opts['save_dir'], 'a.ckpt')
    model = resnet18()
    best_w = torch.load(model_dir)
    model.load_state_dict(best_w['state_dict'])
    if torch.cuda.is_available():
        model = model.cuda()
    val_set = data_loader(dataset, fold=None,is_train=False)
    #val_torchloader =  Data.DataLoader(dataset=val_set, batch_size=1, shuffle=False)

    val_proba = []
    val_true = []
    for feature, label, res_id, pdb_id in val_set:
        # if pdb_id != '7k5l_A':
            # continue
        if pdb_id[-1].islower():
            pdb_id2 = pdb_id + pdb_id[-1]
        else:
            pdb_id2 = pdb_id
        criterion = WeightedMultilabel()
        loss_, proba_ = val_batch(feature, label, model, criterion)
        val_proba.append(proba_.cpu().numpy()[:,1])
        # print(proba_[:,1])
        val_true.append(label.numpy())
        fpr, tpr, thresholds=roc_curve(label.numpy(), proba_.cpu().numpy()[:,1])
        roc_auc=auc(fpr, tpr)
        print(pdb_id, roc_auc)

        #copy useful files into file
        target_folder = os.path.join(dir_opts['predict_dir'], pdb_id2)
        if not os.path.exists(target_folder):
            os.mkdir(target_folder)

        #score.npy
        np.save(os.path.join(target_folder, 'score'), proba_.cpu().numpy()[:,1])
        np.save(os.path.join(target_folder, 'label'), label.cpu().numpy())

        pass

    val_true = np.concatenate(val_true)
    val_proba = np.concatenate(val_proba)
    fpr, tpr, thresholds=roc_curve(val_true, val_proba)
    roc_auc=auc(fpr, tpr)
    print('overall', roc_auc)
    np.save(dir_opts['predict_dir']+'/true_'+output_base, val_true)
    np.save(dir_opts['predict_dir']+'/proba_'+output_base, val_proba)

if __name__ == '__main__':
    dataset_path = '/home/aoli/Documents/RBP_surface_make_database/data/port/case_study.txt'
    # dataset_path = '/home/aoli/Documents/RBP_surface_make_database/data/port/RBP407_09.txt'

    output_base = 'case_study'
    # for fold in range(7,8):
    predict(dataset_path, output_base)
