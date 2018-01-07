import numpy as np
import matplotlib.pyplot as plt

#เก็บ DNA ที่เป็น FASTA เข้ามาใน Dict
x = open("ucsc_hg19.fa","r")
a = 0
data = {}
key = ""
val = ""
for i in x:
    if i[0] == ">" :
        if a!=0 :
            data[key] = val
        key = i[1:].strip()
        val = ""
        a+=1
    else :
        val+=i.strip()
data[key] = val
#เก็บผล RNA Editing
chrList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
           'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
y = open('Human_AG_all_hg19_v2.txt',"r")
rnaEditingPosSense = {}
rnaEditingPosAntisense = {}
rnaEditingPosAll = {}
aluPos = {}
b = 0
for i in y :
    info = i.split()
    key = info[0]
    val = info[1]
    strand = info[3]
    if not key in aluPos :
        aluPos[key] = {}
    if not key in rnaEditingPosSense :
        rnaEditingPosSense[key] = []
        rnaEditingPosAll[key] = []
        rnaEditingPosAntisense[key] = []
    if info[4] == 'intronic' and strand=='+' :
        rnaEditingPosSense[key].append(int(val))
        rnaEditingPosAll[key].append(int(val))
    elif info[4] == 'intronic' and strand=='-' :
        rnaEditingPosAntisense[key].append(int(val))
        rnaEditingPosAll[key].append(int(val))
    aluPos[key][val] = info[6]
rnaEditingPosSense.pop('chrM')
rnaEditingPosSense.pop('chromosome')
rnaEditingPosAntisense.pop('chrM')
rnaEditingPosAntisense.pop('chromosome')
rnaEditingPosAll.pop('chrM')
rnaEditingPosAll.pop('chromosome')
for i in chrList :
    rnaEditingPosSense[i].sort()
    rnaEditingPosAntisense[i].sort()
    rnaEditingPosAll[i].sort()
    
#เตรียมข้อมูล
alldata = {}
alldata['chr19'] = np.zeros((len(data['chr19'])-499,), dtype=np.int)
for i in range(len(rnaEditingPosAll['chr19'])-1) :
    count = 0
    st = rnaEditingPosAll['chr19'][i]
    nx = rnaEditingPosAll['chr19'][i+1]
    delta = nx-st
    if delta < 500 :
        for j in range(st-(500-delta),st+1) :
            alldata['chr19'][j] = 1
    #if i % 100 == 0 : print(i)

X = np.zeros((len(data['chr19'])-499,2),dtype = np.int)
atfirst = 0
for i in range(500) :
    if data['chr19'][i] == "A" or data['chr19'][i] == "a" or data['chr19'][i] == "T" or data['chr19'][i] == "t" :
            atfirst+=1
X[i][1] = atfirst
atdel = False
for i in range(len(X)) :
    X[i][0] = i
    if atdel :
        atfirst-=1
    if data['chr19'][i+499] == "A" or data['chr19'][i+499] == "a" or data['chr19'][i+499] == "T" or data['chr19'][i+499] == "t" :
        atfirst+=1
    X[i][1] = atfirst
    if data['chr19'][i] == "A" or data['chr19'][i] == "a" or data['chr19'][i] == "T" or data['chr19'][i] == "t" :
        atdel = True
    else :
        atdel = False
y = alldata[‘chr19']

#แบ่ง Training Set กับ Test Set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.33, random_state = 0)

#ปรับ Scaling ให้สามารถทำนายได้ง่าย
from sklearn import preprocessing
X_train = preprocessing.scale(X_train)
X_test = preprocessing.scale(X_test)

#สร้าง model แบบ Random Forest Classification
from sklearn.ensemble import RandomForestClassifier
classifier = RandomForestClassifier(n_estimators = 1, criterion = 'entropy', random_state = 0, verbose = 3)
classifier.fit(X_train, y_train)

#ทำนายข้อมูล
y_pred = classifier.predict(X_test)

#สร้าง matrix ดู % ความถูกต้อง
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)

#Plot graph
from matplotlib.colors import ListedColormap
X_set, y_set = X_test, y_test
y_set2 = y_pred
y_set3 = y_diff
X1, X2 = np.meshgrid(np.arange(start = X_set[:, 0].min() - 1, stop = X_set[:, 0].max() + 1),
                     np.arange(start = X_set[:, 1].min() - 1, stop = X_set[:, 1].max() + 1))
plt.xlim(X1.min(), X1.max())
plt.ylim(X2.min(), X2.max())
for i, j in enumerate(np.unique(y_set)):
    plt.scatter(X_set[y_set == j, 0], X_set[y_set == j, 1],
                c = ListedColormap(('red', 'green'))(i), label = j,s=0.005,alpha = 0.3)
plt.title('RNA Editing Test Result [chr19]')
plt.xlabel('Position')
plt.ylabel('ATContent')
plt.legend()
plt.show()
