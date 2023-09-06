import pandas as pd
import numpy as np
import pickle
from sklearn import metrics
from sklearn.metrics import matthews_corrcoef
from sklearn.ensemble import RandomForestClassifier

with open('CNRBind.pkl', 'rb') as f:
    model = pickle.load(f)
    
ori = pd.read_csv(r'TE18.csv')
test_X = ori.iloc[:,:-1]
test_label = ori.iloc[:,-1]
resample_pred = model.predict(np.array(test_X))
y_score = model.predict_proba(np.array(test_X))[:,1]
fpr,tpr,threshold = metrics.roc_curve(test_label, y_score)
roc_auc = metrics.auc(fpr,tpr)
print("The results for CNRBind on TE18")
print('Pre:',metrics.precision_score(test_label,resample_pred))
print("Sn:",metrics.recall_score(test_label,resample_pred))
print('AUC:',roc_auc)
print('MCC:',matthews_corrcoef(test_label,resample_pred))    
print('')

ori = pd.read_csv(r'RB9.csv')
test_X = ori.iloc[:,:-1]
test_label = ori.iloc[:,-1]
resample_pred = model.predict(np.array(test_X))
y_score = model.predict_proba(np.array(test_X))[:,1]
fpr,tpr,threshold = metrics.roc_curve(test_label, y_score)
roc_auc = metrics.auc(fpr,tpr)
print("The results for CNRBind on RB9")
print('Pre:',metrics.precision_score(test_label,resample_pred))
print("Sn:",metrics.recall_score(test_label,resample_pred))
print('AUC:',roc_auc)
print('MCC:',matthews_corrcoef(test_label,resample_pred))    
print('')

ori = pd.read_csv(r'Sub1.csv')
test_X = ori.iloc[:,:-1]
test_label = ori.iloc[:,-1]
resample_pred = model.predict(np.array(test_X))
y_score = model.predict_proba(np.array(test_X))[:,1]
fpr,tpr,threshold = metrics.roc_curve(test_label, y_score)
roc_auc = metrics.auc(fpr,tpr)
print("The results for CNRBind on Sub1")
print('Pre:',metrics.precision_score(test_label,resample_pred))
print("Sn:",metrics.recall_score(test_label,resample_pred))
print('AUC:',roc_auc)
print('MCC:',matthews_corrcoef(test_label,resample_pred))    
print('')

ori = pd.read_csv(r'Sub2.csv')
test_X = ori.iloc[:,:-1]
test_label = ori.iloc[:,-1]
resample_pred = model.predict(np.array(test_X))
y_score = model.predict_proba(np.array(test_X))[:,1]
fpr,tpr,threshold = metrics.roc_curve(test_label, y_score)
roc_auc = metrics.auc(fpr,tpr)
print("The results for CNRBind on Sub2")
print('Pre:',metrics.precision_score(test_label,resample_pred))
print("Sn:",metrics.recall_score(test_label,resample_pred))
print('AUC:',roc_auc)
print('MCC:',matthews_corrcoef(test_label,resample_pred))    
print('')

ori = pd.read_csv(r'Sub3.csv')
test_X = ori.iloc[:,:-1]
test_label = ori.iloc[:,-1]
resample_pred = model.predict(np.array(test_X))
y_score = model.predict_proba(np.array(test_X))[:,1]
fpr,tpr,threshold = metrics.roc_curve(test_label, y_score)
roc_auc = metrics.auc(fpr,tpr)
print("The results for CNRBind on Sub3")
print('Pre:',metrics.precision_score(test_label,resample_pred))
print("Sn:",metrics.recall_score(test_label,resample_pred))
print('AUC:',roc_auc)
print('MCC:',matthews_corrcoef(test_label,resample_pred))    
