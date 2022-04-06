import os
import pandas as pd
import pickle
import re
import numpy as np
from sklearn.metrics import confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt

pos_data = "dataset/pos.fold"
neg_data = "dataset/neg.fold"

print("Welcome to RNAinsecta standalone test run:\n")
def data_prep(path):
 
 print("Reading dataset:\n",path)
 qi = open(path, "r")
 text = qi.read()
 qi.close()

 print("Calculating features:\n")
 matchA = re.findall("[(.][.)(]+[.)]", str(text))
 matchB = re.findall(r"^[AUGC]{10,}", str(text), flags=re.MULTILINE)
 df02 = pd.DataFrame(matchA, matchB)
 df02 = df02.reset_index()
 df02.columns = ["sequence", "fold"]
 os.system('perl trip.pl < ' + path + '> ' + path + 'trip.csv 2>/dev/null')  # miPred processes
 os.system('perl stats.pl <' + path + '> ' + path + 'stats.csv 2>/dev/null') 
 s1 = pd.read_csv(path+"trip.csv", sep='\s+')
 s2 = pd.read_csv(path+"stats.csv", sep='\t')
 s1 = s1.fillna(0)
 s2 = s2.fillna(0)
 s2 = s2.iloc[:, :-1]
 scon = pd.concat([s1, s2], axis=1)
 scon = scon[pd.to_numeric(scon.iloc[:, 0], errors='coerce').notnull()]
 scon = scon.astype(float)
 scon = pd.concat([scon, df02], axis=1)
 scon.rename(columns={"Nmfe": "dG"}, inplace=True)
 scon["stem"] = scon["fold"].str.extract('(\(.*\))')
 scon["nstem"] = scon["stem"].apply(len)
 scon["stemte"] = scon["fold"].str.findall('(\(\.+\))')
 scon["nstemte"] = scon["stemte"].apply(len)
 scon["stemte"] = scon["fold"].str.extract('(\(\.+\))')
 scon["lenstemte"] = scon["stemte"].apply(len)

 scon["MFE1"] = scon["dG"] / scon["%G+C"]
 scon["stem0"] = scon["fold"].str.findall('(\({3,})')
 scon["n_stems"] = scon["stem0"].apply(len)
 scon["stem1"] = scon["fold"].str.findall('(\()')
 scon["total_base"] = scon["stem1"].apply(len)
 scon["MFE2"] = scon["dG"] / scon["n_stems"]
 scon["MFE3"] = scon["dG"] / scon["nstemte"]
 scon["MFE4"] = scon["dG"] / scon["total_base"]
 scon["avg_bp"] = scon["total_base"] / scon["n_stems"]

 scon1 = scon[
  ["A...", "A..(", "A.(.", "A.((", "A(..", "A(.(", "A((.", "A(((", "G...", "G..(", "G.(.", "G.((", "G(..",
   "G(.(",
   "G((.", "G(((", "C...", "C..(", "C.(.", "C.((", "C(..", "C(.(", "C((.", "C(((", "U...", "U..(", "U.(.",
   "U.((",
   "U(..", "U(.(", "U((.", "U(((", "Len", "A", "C", "G", "U", "G+C", "A+U", "AA", "AC", "AG", "AU", "CA",
   "CC",
   "CG",
   "CU", "GA", "GC", "GG", "GU", "UA", "UC", "UG", "UU", "%A", "%C", "%G", "%U", "%G+C", "%A+U", "%AA", "%AC",
   "%AG",
   "%AU", "%CA", "%CC", "%CG", "%CU", "%GA", "%GC", "%GG", "%GU", "%UA", "%UC", "%UG", "%UU", "pb", "Npb",
   "mfe",
   "dG", "Q", "NQ", "D", "ND", "nstem", "MFE1", "MFE2", "MFE3", "MFE4", "total_base", "n_stems", "avg_bp"]]

 b = scon[["sequence", "Len", "mfe", "fold", "avg_bp"]]
 scon1 = scon1.abs()
 return scon1


pos = data_prep(pos_data)
neg = data_prep(neg_data)

pos.loc[:, 'SEQ'] = int(1)
neg.loc[:, 'SEQ'] = int(0)

tot = pd.concat([pos,neg])
tot = tot.abs()

X = tot.iloc[:, :-1].values
y = tot.iloc[:, -1].values

print("Running predictive model:\n")

svm_path = "ML_models/SVM.pkl"
rf_path = "ML_models/RF.pkl"

def testing(model):
 classifier = pickle.load(open(model, 'rb'))
 y_pred = classifier.predict(X)
 y_pred_proba = classifier.predict_proba(X)
 y_proba = y_pred_proba[:, 1]
 print("Running complete for\n",model)
 return y_pred, y_proba;

def eval(y_pred):
 coX = confusion_matrix(y, y_pred)
 neg = coX[0]
 pos = coX[1]
 tot = np.append(neg, pos)
 tot = np.reshape(tot, (1, 4))
 df = pd.DataFrame(tot, columns=["TN", "FP", "FN", "TP"])
 df["acc"] = (df["TP"] + df["TN"]) / (df["TP"] + df["FP"] + df["TN"] + df["FN"])
 df["SP"] = df["TN"] / (df["TN"] + df["FP"])
 df["SN"] = df["TP"] / (df["TP"] + df["FN"])
 df["MCC"] = ((df["TP"] * df["TN"]) + (df["FP"] * df["FN"])) / np.sqrt(
  (df["TN"] + df["FP"]) * (df["TN"] + df["FN"]) * (df["TP"] + df["FP"]) * (df["TP"] + df["FN"]))
 df["p"] = df["TP"] / (df["TP"] + df["FP"])
 df["F1"] = 2 * (df["SN"] * df["p"]) / (df["SN"] + df["p"])
 return df

svm_pred, svm_pred_proba = testing(svm_path)
svm_res = eval(svm_pred)

rf_pred, rf_pred_proba = testing(rf_path)
rf_res = eval(rf_pred)

res = pd.concat([rf_res, svm_res], ignore_index=True)
res = res.rename(index={0: 'RF', 1:'SVM'})

os.system("mkdir results:")
res.to_csv("results/results.csv")
print("Prediction Results:\n",res)
print("\n\nThe files are stored in the results directory.")
print("\nThank you for using RNAinsecta. Kindly send your queries or issues to adhiraj@iitg.ac.in")
#AUC_ROC
svm_fpr, svm_tpr, threshold = roc_curve(y, svm_pred_proba)
svm_auc = auc(svm_fpr, svm_tpr)

rf_fpr, rf_tpr, threshold = roc_curve(y, rf_pred_proba)
rf_auc = auc(rf_fpr, rf_tpr)

plt.figure(figsize=(5, 5), dpi=800)
plt.plot(svm_fpr, svm_tpr, linestyle='-', label='SVM (auc = %0.3f)' % svm_auc)
plt.plot(rf_fpr, rf_tpr, linestyle='-', label='RF (auc = %0.3f)' % rf_auc)
plt.xlabel('False Positive Rate -->')
plt.ylabel('True Positive Rate -->')
plt.legend()
plt.show()
fig = plt.gcf()
fig.set_size_inches((10, 6), forward=False)
fig.savefig('results/ROC.png', format='png', dpi=800)