

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import savefig
from scipy.io import arff
import ntpath
import glob
import os
import math
from sklearn.model_selection import train_test_split
from random import randrange
from itertools import combinations as comb
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.metrics import accuracy_score
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import Lasso
from sklearn.model_selection import cross_val_score
from sklearn.metrics import accuracy_score, f1_score
from sklearn.metrics import classification_report
from svglib.svglib import svg2rlg
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import scikitplot as skplt
import matplotlib.pyplot as plt
import sklearn.metrics as mt
from google.colab import drive
drive.mount('/content/gdrive')

gdrivePath = "gdrive" + os.sep +"My Drive"
resultsPath = gdrivePath + os.sep + "Stage Classification" + os.sep + "Figures"

#read the dataset and the labels
df = pd.read_csv(gdrivePath +os.sep+"Stage Classification"+os.sep+"Data"+os.sep+"stage_vs_others_dge.csv", index_col=0)
df_label = pd.read_csv(gdrivePath +os.sep+"Stage Classification"+os.sep+"Data"+os.sep+"integrated_pd.csv", index_col=0)
df_label = df_label['stageall']
df_label = pd.DataFrame(df_label)

df.shape

df_label.value_counts()

#equalization number of samples in each class
from imblearn.over_sampling import SMOTE
smote = SMOTE(sampling_strategy='all')
X_sm, Y_sm = smote.fit_resample(df, df_label)
Y_sm.value_counts()

#Split the data
X_train, X_test, y_train, y_test = train_test_split(X_sm, Y_sm, test_size=0.2, random_state=0, stratify = Y_sm)


names = ["Quadratic_Discriminant_Analysis","Extra_Trees","Linear_SVM",  "Random_Forest",  "3_Nearest_Neighbors"    ]

def classify(X_train, X_test, y_train, y_test):
 

  names = ["Quadratic_Discriminant_Analysis","Extra_Trees","Linear_SVM",  "Random_Forest",  "3_Nearest_Neighbors"    ]        

  classifiers = [
    QuadraticDiscriminantAnalysis(),
    ExtraTreesClassifier(random_state=0, criterion = 'entropy'),  
    SVC(kernel="linear", C=0.025, probability=True),
    RandomForestClassifier(random_state=0, criterion = 'entropy'),
    KNeighborsClassifier(3)
    ]

  acc_scores = []
  train_scores = []
  val_scores = []
  test_predictions = []
  for name, clf in zip(names, classifiers):
      clf.fit(X_train, y_train)
      acc_score = clf.score(X_test, y_test)
      train_score = clf.score(X_train, y_train)
      val_score = cross_val_score(clf, X_train, y_train, cv=5).mean()
      acc_scores.append(acc_score)
      train_scores.append(train_score)
      val_scores.append(val_score)
      test_predictions.append(clf.predict(X_test))

  dfs = pd.DataFrame()
  dfs['Method'] = names
  dfs['Train '] = train_scores
  dfs['5-Fold Validation'] = val_scores
  dfs['Test'] = acc_scores
  
  return dfs, test_predictions


dfs,pr = classify(X_train, X_test, y_train, y_test)


for i in range(0,5):
  print("\nMCC Score \n",names[i], mt.matthews_corrcoef(y_test, predictions[i]))

for i in range(0,5):
print("\nROC AUC Score\n", mt.roc_auc_score(y_test, predictions[i],labels=[1,2,3,4], average='macro',multi_class='ovo'))

for i in range(0,5):
  print("\nPrecision Score \n",names[i], mt.precision_score(y_test, predictions[i],labels=[1,2,3,4], average='macro'))

for i in range(0,5):
  print("\nRecall Score \n",names[i], mt.recall_score(y_test, predictions[i],labels=[1,2,3,4], average='macro'))


LabelDict = { 1 : "Stage-1", 2: "Stage-2", 3:"Stage-3", 4:"Stage-4"}
cm_mlp = confusion_matrix(y_test, y_pred)
disp = ConfusionMatrixDisplay(confusion_matrix=cm_mlp, display_labels=  [LabelDict[cat] for cat in clf.classes_])
disp.plot()


target_names = ['Stage-1', 'Stage-2', 'Stage-3','Stage-4' ]

for i in range(0,5):
  print('\033[91m' '\033[1m' + dfss["Yöntem"][i] + '\033[0m')
  print((classification_report(y_test, predictions[i], target_names=target_names)))




names = ["QDA","Extra Trees","SVM",  "Random Forest",  "3 Nearest Neighbors"    ] 
scores = [79, 77, 72, 75, 61]
plt.figure(figsize=(15, 8))
plt.barh(names, scores)
plt.ylabel("model")
plt.xlabel("accuracy (%)")
plt.title("Test Accuracy Results")
image_name = resultsPath + os.sep +'test_acc_bar_svg'
plt.savefig(image_name, format='svg', dpi = 100)



names = ["QDA","Extra Trees","SVM",  "RF",  "3-NN" ] 
scores = [79, 77, 72, 75, 61]
without_smote = [31,52,43,49,43]
  
X_axis = np.arange(len(names))
  
plt.bar(X_axis - 0.2, scores, 0.4, label = 'with SMOTE')
plt.bar(X_axis + 0.2, without_smote, 0.4, label = 'without SMOTE')
  
plt.xticks(X_axis, names)
plt.xlabel("Models")
plt.ylabel("Test Accuracy (%)")
plt.title("SMOTE Effect")
plt.legend()
plt.savefig(image_name, format='png', dpi=1200)


