import os, sys
import torch
import numpy as np
import scipy.io as scio
from torch import optim
import torch.nn.functional as F

from nn_module.resnet18 import resnet18
from default_config.config import config
from default_config.dir_options import dir_opts
from nn_module.function import f1, auc, save_ckpt, WeightedMultilabel,adjust_learning_rate
from data_loader.data_loader import data_loader
import torch.utils.data as Data

os.environ["CUDA_VISIBLE_DEVICES"] = "1"
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print('Model is working on device', device)

def train_batch(protein_data, site_label, model, criterion, optimizer):
    model.train()
    if torch.cuda.is_available():
        protein_data = protein_data.cuda()
        site_label = site_label.cuda()
    output = model(protein_data)
    loss = criterion(output,site_label)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    loss = loss.data.item()
    return loss, F.softmax(output, dim=1)

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

def train(dataset, fold):
    learning_rate = config['lr']
    num_epoches = config['epochs']
    train_set = data_loader(dataset, fold, is_train=True)
    train_set.stastic()
    train_torchloader = Data.DataLoader(dataset=train_set, batch_size=config['batch_size'], shuffle=True)

    val_set = data_loader(dataset, fold, is_train=False)
    val_set.stastic()
    val_torchloader =  Data.DataLoader(dataset=val_set, batch_size=config['batch_size'], shuffle=False)


    model = resnet18()
    if torch.cuda.is_available():
        model = model.cuda()
    criterion = WeightedMultilabel()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    if not os.path.exists(dir_opts['save_dir']):
        os.mkdir(dir_opts['save_dir'])
    #train model
    loss_on_epoch = []; val_loss_on_epoch = []; auc_on_epoch = []; val_auc_on_epoch = []
    auc_tmp = 0
    for epo in range(num_epoches):

        loss_on_batches = []; proba = []; true = []
        batch_number = 0
        for feature, label in train_torchloader:
            loss_, proba_ = train_batch(feature, label, model, criterion, optimizer)
            loss_on_batches.append(loss_)
            proba.append(proba_.cpu())
            true.append(label)
            if batch_number % 100 == 99:
                out_str = '\rEpoch:{:0>2d}, loss:{:.3f}'.format(epo, np.mean(np.array(loss_on_batches)[-100:]))
                sys.stdout.write(out_str)
                sys.stdout.flush()
            batch_number += 1

        #evaluate on training dataset after one epoch
        train_loss_on_batches = []; train_proba = []; train_true = []
        for feature, label in train_torchloader:
            loss_, proba_ = val_batch(feature, label, model, criterion)
            train_loss_on_batches.append(loss_)
            train_proba.append(proba_.cpu())
            train_true.append(label)

        train_loss_on_batches = np.average(np.array(train_loss_on_batches))
        train_proba = torch.vstack(train_proba)
        train_true = torch.hstack(train_true)
        auc_train = auc(train_true, train_proba)
        loss_on_epoch.append(train_loss_on_batches)
        auc_on_epoch.append(auc_train)

        #evaluate on validation dataset after one epoch
        val_loss_on_batches = []; val_proba = []; val_true = []
        for feature, label in val_torchloader:
            loss_, proba_ = val_batch(feature, label, model, criterion)
            val_loss_on_batches.append(loss_)
            val_loss_on_epoch.append(loss_)
            val_proba.append(proba_.cpu())
            val_true.append(label)

        val_loss_on_batches = np.average(np.array(val_loss_on_batches))
        val_proba = torch.vstack(val_proba)
        val_true = torch.hstack(val_true)
        auc_val = auc(val_true, val_proba)
        val_loss_on_epoch.append(val_loss_on_batches)
        val_auc_on_epoch.append(auc_val)


        out_str = '\rEpoch:{:0>2d}, train_loss:{:.3f}, train_auc:{:.3f}, val_loss:{:.3f}, val_auc:{:.3f}'.format(epo,
                                                                              train_loss_on_batches, auc_train, val_loss_on_batches, auc_val)
        sys.stdout.write(out_str+'\n')
        sys.stdout.flush()

        #保存模型
        if auc_val > auc_tmp:
            auc_tmp = auc_val
            state = {"state_dict": model.state_dict(), "epoch": epo}
            save_ckpt(state,fold)

        if epo in config['stage_epoch']:
            learning_rate /= config['lr_decay']
            adjust_learning_rate(optimizer, learning_rate)
            model.load_state_dict(torch.load(os.path.join(dir_opts['save_dir'], 'best_'+ str(fold)+'.ckpt'))['state_dict'])
            print('change learning rate ')
    #保存曲线
    loss_on_epoch = np.array(loss_on_epoch)
    val_loss_epoch = np.array(val_loss_on_epoch)
    auc_on_epoch = np.array(auc_on_epoch)
    val_auc_on_epoch = np.array(val_auc_on_epoch)

    scio.savemat(os.path.join(dir_opts['save_dir'],'history'+str(fold)+'.mat'), {'train_loss': loss_on_epoch,'train_auc': auc_on_epoch,
                                                                    'val_loss': val_loss_epoch,'val_auc':val_auc_on_epoch})

def predict(dataset, fold):
    model = resnet18()
    best_w = torch.load(os.path.join(dir_opts['save_dir'], 'best_'+ str(fold)+'.ckpt'))
    model.load_state_dict(best_w['state_dict'])
    if torch.cuda.is_available():
        model = model.cuda()
    val_set = data_loader(dataset, fold, is_train=False)
    val_torchloader =  Data.DataLoader(dataset=val_set, batch_size=config['batch_size'], shuffle=False)
    criterion=WeightedMultilabel()

    val_proba = []
    val_true = []
    for feature, label in val_torchloader:
        loss_, proba_ = val_batch(feature, label, model, criterion)
        val_proba.append(proba_.cpu().numpy()[:,1])
        val_true.append(label.numpy())
    val_true = np.concatenate(val_true)
    val_proba = np.concatenate(val_proba)

    np.save(dir_opts['cross_dir']+'/true_fold'+str(fold), val_true)
    np.save(dir_opts['cross_dir']+'/proba_fold'+str(fold), val_proba)

if __name__ == '__main__':
    dataset = dir_opts['dataset']
    for fold in range(0, 10):
        print('**************', fold)
        train(dataset, fold)