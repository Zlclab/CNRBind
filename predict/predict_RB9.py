import copy
import os
import networkx as nx
import pickle
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from sklearn import metrics
from sklearn.metrics import matthews_corrcoef
from sklearn.ensemble import RandomForestClassifier

name = []#########A list for all RNA name
sequence = []#########A list for all RNA sequence
with open('..\\data_cashe\\RB9.txt', 'r') as f:
    lines = f.readlines()
for i in range(len(lines)):
    if i % 3 == 0:
        name.append(lines[i][1:-1])
    if i % 3 == 1:
        sequence.append(lines[i][:-1])

######################ASA was obtained by RNAsol
ori = pd.read_excel(r'RB9_ASA&label.xlsx')
ASA = ori.iloc[:,0]
label = ori.iloc[:,1]

seq_len = []#########A list for all RNA length
for i in range(len(sequence)):
    seq_len.append(len(sequence[i]))

def Kmers_funct(seq, size): 
    return [seq[x:x+size] for x in range(len(seq) - size + 1)]

one_mers=[]
for i in range(len(sequence)):
    myseq = sequence[i]
    one_mers.append(Kmers_funct(myseq, size=1))

all_AA=[]#########A long list for all amino acids
for i in range(len(one_mers)):
    all_AA.extend(one_mers[i])

Name = []
for i in range(len(sequence)):
    for j in range(len(sequence[i])):
        Name.append(name[i])
        
########################Start building nucleotide information-based feature vector

########################site location coefficient (SL)
sl = []
for i in range(len(sequence)):
    one = Kmers_funct(sequence[i],1)
    for j in range(len(one)):
        sl.append((j+1)/len(one))
        
########################nucleotide position-specific coefficient (NPS)     
nps = []
for i in range(len(sequence)):
    if len(sequence[i])%2 == 1:
        for j in range(int((len(sequence[i])+1)/2),0,-1):
            nps.append(j/((len(sequence[i])+1)/2))
        for j in range(2,int((len(sequence[i])+3)/2)):
            nps.append(j/((len(sequence[i])+1)/2))
    else:
        for j in range(int(len(sequence[i])/2),0,-1):
            nps.append(j/(len(sequence[i])/2))
        for j in range(1,int((len(sequence[i])+2)/2)):
            nps.append(j/(len(sequence[i])/2))
            
########################density (DST)
dst=[]
for i in range(len(sequence)):
    for j in range(len(sequence[i])):
        a=sequence[i][:j+1]
        dst.append(a.count(a[j])/(j+1))

########################frequency of occurrence (Fre)
fre = []
for i in range(len(sequence)):
    one = Kmers_funct(sequence[i],1)
    for j in range(len(one)):
        fre.append(one.count(one[j])/len(one))

################EIIP+CB+MM+pKa+RFHC+CS+BE+BER (See the RNA.xlsx file for details.)
ori1 = pd.read_excel(r'RNA.xlsx')
nuc = ori1.iloc[:,0]
inf = ori1.iloc[:,1:]
info = ori1.set_index('nuc').T.to_dict('list')

iii=[]
for i in range(len(all_AA)):
    iii.append(info.get(all_AA[i]))
df1 = pd.DataFrame(iii)

############################NDC and NDS
path = "..\\pdb\\"###### location of PDB files
df_empty = pd.DataFrame()
NDS = []###Euclidean distance sum 
NDC = []###Euclidean distance to geometric center
COSS = []###cosine distance  sum 
COSC = []###cosine distance  to geometric center
CHEBS = []###Chebyshev distance sum 
CHEBC = []###Chebyshev distance to geometric center
for n in range(len(name)):    
    x1 = []
    y1 = []
    z1 = []
    for line in open(path + name[n] +'.pdb'):
        list = line.split()
        if list[0] == 'ATOM'or list[0] == 'HETATM':
            if list[2] == "C1'":
                x1.append(float(list[6]))
                y1.append(float(list[7]))
                z1.append(float(list[8]))
    zhixin = np.array([sum(x1)/len(x1),sum(y1)/len(y1),sum(z1)/len(z1)])
    
    for i in range(len(x1)):
        nds = []
        coss = []
        chebs = []
        for j in range(len(x1)):
            if i!=j:
                a = np.array([x1[i],y1[i],z1[i]])
                b = np.array([x1[j],y1[j],z1[j]])
                X = np.vstack([a,b])
                c = pdist(X)
                nds.append(c[0])
                coss.append(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))
                chebs.append(np.abs(a-b).max())
        NDS.append(sum(nds)/len(x1))
        COSS.append(sum(coss)/len(x1))
        CHEBS.append(sum(chebs)/len(x1))
        d = np.array([x1[i],y1[i],z1[i]])
        Y = np.vstack([d,zhixin])
        e = pdist(Y)
        NDC.append(e[0])
        COSC.append(np.dot(d,zhixin)/(np.linalg.norm(d)*np.linalg.norm(zhixin)))
        CHEBC.append(np.abs(d-zhixin).max())

########################################Summarize
nuc_inf_vec = {'ASA':ASA,'sl':sl,'dst':dst,'fre':fre,'nps':nps,'BE1':df1.iloc[:,11],'BE2':df1.iloc[:,12],
               'BE3':df1.iloc[:,13],'BE4':df1.iloc[:,14],'BER1':df1.iloc[:,15],'BER2':df1.iloc[:,16],
               'BER3':df1.iloc[:,17],'BER4':df1.iloc[:,18],'RFHC1':df1.iloc[:,4],'RFHC2':df1.iloc[:,5],
               'RFHC3':df1.iloc[:,6],'CS1':df1.iloc[:,7],'CS2':df1.iloc[:,8],'CS3':df1.iloc[:,9],
               'CS4':df1.iloc[:,10],'CB':df1.iloc[:,1],'EIIP':df1.iloc[:,0],'pKa':df1.iloc[:,3],
               'MM':df1.iloc[:,2],'NDS':NDS,'NDC':NDC,'COSS':COSS,'COSC':COSC,'CHEBS':CHEBS,'CHEBC':CHEBC,
               }
df_nuc_inf = pd.DataFrame(nuc_inf_vec)
df_nuc_inf

######################################Construct complex networks
df_empty = pd.DataFrame()

for n in range(len(name)):    
    x = []
    y = []
    z = []
    NO_aminoacid = []
    for line in open(path + name[n] +'.pdb'):
        list = line.split()
        if list[0] == 'ATOM'or list[0] == 'HETATM':
            NO_aminoacid.append(list[5])
            x.append(float(list[6]))
            y.append(float(list[7]))
            z.append(float(list[8]))
    number_of_atom = len(x)
    number_of_aminoacid = 1
    Rev_NO_aminoacid = [1]*number_of_atom
    for i in range(2,len(x)):
        if NO_aminoacid[i]!=NO_aminoacid[i-1]:
            number_of_aminoacid = number_of_aminoacid +1
            for j in range(i,len(x)):
                Rev_NO_aminoacid[j] = number_of_aminoacid
            
    contact = np.zeros((number_of_aminoacid, number_of_aminoacid)).astype('int64')
    for i in range(len(x)):
        for j in range(len(x)):
            if abs(Rev_NO_aminoacid[i]-Rev_NO_aminoacid[j]) > 1:
                a=[x[i],y[i],z[i]]
                b=[x[j],y[j],z[j]]
                X=np.vstack([a,b]) 
                d_ij = pdist(X)
                if d_ij <= 8 :
                    contact[Rev_NO_aminoacid[i]-1][Rev_NO_aminoacid[j]-1] = 1
    G=nx.Graph(contact)
###################################calculate DG
    degree = []
    degree.extend(nx.degree_centrality(G).values())

###################################calculate CL
    closeness = []
    closeness.extend(nx.closeness_centrality(G).values())

###################################calculate BC
    betweenness = []
    betweenness.extend(nx.betweenness_centrality(G).values())

###################################Find the adjacent nodes of each site
    son_of_site = []
    for i in range(len(contact)):
        site = []
        for j in range(len(contact)):
            if contact[i][j] == 1:
                site.append(j+1)
        son_of_site.append(site)
    
#calculate DG & CL & BC for adjacent nodes 
    a1 = []
    a2 = []
    a3 = []
    for i in range(len(son_of_site)):
        b1 = []
        b2 = []
        b3 = []
        for j in range(len(son_of_site[i])):
            b1.append(degree[son_of_site[i][j]-1])
            b2.append(closeness[son_of_site[i][j]-1])
            b3.append(betweenness[son_of_site[i][j]-1])
        a1.append(b1)
        a2.append(b2)
        a3.append(b3)
       
    
    degree_of_son = []
    closeness_of_son = []
    betweenness_of_son = []
    for i in range(len(a1)):
        degree_of_son.append(sum(a1[i]))
        closeness_of_son.append(sum(a2[i]))
        betweenness_of_son.append(sum(a3[i]))
    info_ = []
    for j in range(len(Name)):
        if Name[j] == name[n]: 
            info_.append(df_nuc_inf.iloc[j,:])
    info = pd.DataFrame(info_)
    info.reset_index(drop = True,inplace = True)
    
    info2 = []
    for i in range(len(son_of_site)):
        info1 = []
        for j in range(len(son_of_site[i])):
            info1.append(info.loc[son_of_site[i][j]-1])
        info2.append(info1)
    
    info51 = []
    info61 = []
    info52 = []
    for i in range(len(info2)):
        info31 = []
        info41 = []
        info32 = []
        for j in range(len(info2[i])):
            for k in range(len(info2[i][j])):
                info31.append(info2[i][j][k] * a1[i][j])
                info41.append(info2[i][j][k] * a2[i][j])
                info32.append(info2[i][j][k] * a3[i][j])
        info51.append(info31)
        info61.append(info41)
        info52.append(info32)
    
    info7 = []
    for i in range(len(info51)):
        x1 = np.array(info51[i])
        x2 = int(len(x1)/30)
        x3 = np.reshape(x1, (x2,30))
        info7.append(np.sum(x3, axis=0))

    info8 = []
    for i in range(len(info61)):
        y1 = np.array(info61[i])
        y2 = int(len(y1)/30)
        y3 = np.reshape(y1, (y2,30))
        info8.append(np.sum(y3, axis=0))
        
    info9 = []
    for i in range(len(info52)):
        x1 = np.array(info52[i])
        x2 = int(len(x1)/30)
        x3 = np.reshape(x1, (x2,30))
        info9.append(np.sum(x3, axis=0))
   
    e=np.hstack((info7,info8,info9))
    df_e = pd.DataFrame(e)
    df_empty = pd.concat([df_empty,df_e],ignore_index=True)

########################################normalized to the range of [0,1]
seq_len_sum = []
for i in range(len(seq_len)):
    if i == 0:
        seq_len_sum.append(seq_len[i])
    else:
        seq_len_sum.append(sum(seq_len[:i+1]))
seq_len_sum.insert(0,0)

for i in range(len(seq_len)):
    df = df_empty.iloc[seq_len_sum[i]:seq_len_sum[i+1],:]
    for j in df.columns:
        Max = np.max(df[j])
        Min = np.min(df[j])
        df[j] = (df[j] - Min)/(Max - Min)
    df_empty.iloc[seq_len_sum[i]:seq_len_sum[i+1],:] = df

#######################################sliding window
def yilie(array):
    a = []
    for i in range(len(array)):
        for j in range(len(array[0])):
            a.append(array[i][j])
    return a

def window(data,w):
    data_up = data
    for i in range(1,11):
        data_up = np.row_stack((data[i],data_up))
    
    data_down = data_up
    for i in range(len(data)-2,len(data)-12,-1):
         data_down = np.row_stack((data_down,data[i]))

    win = []
    for i in range(len(data)):
        win.append(yilie(data_down[10+i-w:11+i+w,:]))
        
    return pd.DataFrame(win)


def stack(ww):
    empty = pd.DataFrame(data=None,columns=range((2*ww+1)*len(df_empty.iloc[0])),index=range(len(df_empty)))
    for i in range(len(seq_len)):
        df = np.array(df_empty.iloc[seq_len_sum[i]:seq_len_sum[i+1],:])
        empty.iloc[seq_len_sum[i]:seq_len_sum[i+1],:] = window(df,ww)
    return empty

########################################predict

with open('CNRBind.pkl', 'rb') as f:
    model = pickle.load(f)
    
test_X = stack(5)
test_label = label
resample_pred = model.predict(np.array(test_X))
y_score = model.predict_proba(np.array(test_X))[:,1]
fpr,tpr,threshold = metrics.roc_curve(test_label, y_score)
roc_auc = metrics.auc(fpr,tpr)
print('Pre:',metrics.precision_score(test_label,resample_pred))
print("Sn:",metrics.recall_score(test_label,resample_pred))
print('AUC:',roc_auc)
print('MCC:',matthews_corrcoef(test_label,resample_pred))    
