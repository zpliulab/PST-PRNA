import os
from default_config.dir_options import dir_opts
import torch
from torch import nn
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score

def save_ckpt(state, fold):
    if not os.path.exists(dir_opts['save_dir']):
        os.mkdir(dir_opts['save_dir'])
    current_w = os.path.join(dir_opts['save_dir'], 'best_'+ str(fold)+'.ckpt')
    torch.save(state, current_w)

#       true   0  1
#predict  0    TN  FP
#         1    FN  TP
def precision(true, predict):
    #sen = TP / (FN+TP)
    true = true.to(torch.float32)
    predict = torch.softmax(predict,dim=1)
    predict = predict[:,1]
    tmp = torch.mul(true,predict)
    TP = torch.sum(tmp.round())
    predict_P = torch.sum(predict.round())
    pre = TP / (predict_P + 0.00001)
    return pre

def sensitivity(true, predict):
    #pre = TP / (FP+TP)
    true = true.to(torch.float32)
    predict = torch.softmax(predict,dim=1)
    predict = predict[:,1]
    tmp = torch.mul(true,predict)
    TP = torch.sum(tmp.round())
    P_num = torch.sum(true)
    sen = TP / (P_num + 0.00001)
    return sen

def specificity(true, predict):
    true = 1-true.to(torch.float32)
    predict = torch.softmax(predict,dim=1)
    predict = predict[:,0]
    tmp = torch.mul(true,predict)
    TN = torch.sum(tmp.round())
    N_num = torch.sum(true)
    spe = TN / (N_num + 0.00001)
    return spe

def f1(true, predict):
    #f1 = 2*pre*recall / pre+recall
    predict = F.softmax(predict,dim=1)
    sen = sensitivity(true, predict)
    pre = specificity(true, predict)

    f_measure = 2*pre*sen / (sen + pre + 0.00001)
    return f_measure

def auc(true, predict):
    predict = predict.detach().numpy()
    true = true.numpy()
    return roc_auc_score(true,predict[:,1])

class WeightedMultilabel(nn.Module):
    def __init__(self, weights=None):
        super(WeightedMultilabel, self).__init__()
        self.weights = weights

        if self.weights == None:
            self.cerition = nn.CrossEntropyLoss()
        else:
            self.cerition = nn.CrossEntropyLoss(reduction='none')
            self.weights = weights

    def forward(self, predict, true):
        if self.weights == None:
            loss = self.cerition(predict, true)
            return loss
        else:
            loss = self.cerition(predict, true)
            loss = loss * (true*self.weights+1)
            return torch.mean(loss)

def adjust_learning_rate(optimizer, lr):
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr
    return lr
