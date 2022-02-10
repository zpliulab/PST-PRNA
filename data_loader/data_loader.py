from default_config.dir_options import dir_opts
import os, random
import numpy as np
import torch
import torchvision.transforms as transform
from sklearn.model_selection import KFold
import _pickle as cPickle

class data_loader():
    def __init__(self, dataset, fold, is_train=True):
        self.is_train = is_train
        self.protein_list = []
        with open(dataset, 'r') as pid:
            for line in pid.readlines():
                self.protein_list.append(line.strip())
        random.seed(2021)
        random.shuffle(self.protein_list)
        kf = KFold(n_splits=2)
        if fold!= None:
            for index, (train_index, test_index) in enumerate(kf.split(self.protein_list)):
                if index == fold:
                    if is_train:
                        self.protein_list = [self.protein_list[i] for i in train_index]
                    else:
                        self.protein_list = [self.protein_list[i] for i in test_index]
        self.data_all = []
        self.label_all = []
        success = 0
        for pdb in self.protein_list:
            try:
                with open(os.path.join(dir_opts['data_label'], pdb), 'rb') as pid:
                    data0 = cPickle.load(pid)
                for index in range(data0['data'].shape[0]):
                    data = data0['data'][index]
                    label = data0['label'][index]
                    data = data.transpose([2,0,1])
                    self.data_all.append(torch.tensor(data, dtype=torch.float32))
                    self.label_all.append(torch.tensor(label, dtype=torch.long))
                    success += 1
            except:
                print(pdb,'not loaded')

        if is_train:
            print('training dataset:', success)
        else:
            print('validation dataset:', success)

        self.transform = transform.Compose([transform.RandomRotation(180),])

    def __getitem__(self, item):
        data = self.data_all[item]
        if self.is_train:
            data = self.transform(data)

        label = self.label_all[item]
        return data, label

    def __len__(self):
        return len(self.data_all)

    def stastic(self):
        print('positive samples:',np.sum(self.label_all))
        print('negative samples', len(self.label_all)- np.sum(self.label_all))