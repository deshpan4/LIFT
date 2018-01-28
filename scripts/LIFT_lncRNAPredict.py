import argparse
import sys, getopt
from optparse import OptionParser
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import os
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2 import robjects as ro
import sys, getopt
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import Lasso
from rpy2.robjects.vectors import ListVector

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-tr', '--training', help="input training set feature matrix")
    parser.add_argument('-te', '--test', help="input test set feature matrix")
    parser.add_argument('-o', '--output', help="output test set filename")
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

pandas2ri.activate()
ro.conversion.py2ri = ro.numpy2ri
ro.numpy2ri.activate()
df1 = pd.read_csv(args.training)
df2 = pd.read_csv(args.training)
df3 = pd.read_csv(args.test)
df31 = np.asarray(pd.read_csv(args.test))
Xtrain = np.asarray(df1.iloc[:, :-1])
Xtest = np.asarray(df3.iloc[:, :-1])
ytr=df2.iloc[:,-1].values
yte=np.random.randint(1,3,len(df3))
irf = importr("iRF")
auc = importr("AUC")
base = importr("base")
graphics = importr("graphics")
grdevices = importr("grDevices")
R = ro.r
Xtr_mat1 = numpy2ri(Xtrain)
Xte_mat1 = numpy2ri(Xtest)
ytr_mat = ro.FactorVector(ytr)
yte_mat = ro.FactorVector(yte)
Xtr_mat = r.assign("bar", Xtr_mat1)
Xte_mat = r.assign("bar", Xte_mat1)
tempyte_mat = ytr_mat
ncol=robjects.r('ncol')
rep=robjects.r('rep')
p1=ncol(df2)
p=p1[0]
selprob = rep(1/p,p)
rf = robjects.r('list()')
b = irf.randomForest(Xtr_mat, ytr_mat, Xte_mat, base.sample(yte_mat), selprob, ntree=400)
print('Prediction finished')