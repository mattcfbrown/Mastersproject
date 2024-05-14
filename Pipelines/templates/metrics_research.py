#This script will take in individual data and return a metric (NOT A TEXT FILE)

import numpy as np
import sklearn.metrics as metrics
import sys

#Firstly import the data needed
given = np.genfromtxt(sys.argv[1], delimiter = ",")
original = np.genfromtxt(sys.argv[2], delimiter = ",")
original = [x for row in original for x in row]
given = [x for row in given for x in row]

if max(given) == 2:
    Matthew = None
    AUROC = None
    print(Matthew)
    print(AUROC)
else:
    Matthew = metrics.matthews_corrcoef(original,given)
    fpr, tpr, threshold = metrics.roc_curve(original,given)
    AUROC = metrics.auc(fpr,tpr)
    print(Matthew)
    print(AUROC)