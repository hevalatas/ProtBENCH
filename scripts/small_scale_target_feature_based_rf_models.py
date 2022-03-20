# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from datetime import datetime

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, precision_score, recall_score, accuracy_score, matthews_corrcoef


def RF_model(tr,ts,prot_ft,chembl_id,fold):
    start = datetime.now()
    
    df_prot_ft = pd.read_csv(r"..\datasets\small_scale\feature_vectors\{}.tsv".format(prot_ft),sep="\t")

    train_ft = tr.merge(df_prot_ft,on="target_id")
    X_train = train_ft.iloc[:,2:]
    y_train = train_ft["bioactivity"]
    
    test_ft = ts.merge(df_prot_ft,on="target_id")
    X_test = test_ft.iloc[:,2:]
    y_test = test_ft["bioactivity"]
       
    #train model
    clf = RandomForestClassifier(n_estimators=200,max_features="sqrt",random_state=42)
    clf.fit(X_train, y_train)        
           
    # #save model
    # with open(r"..\models\small_scale\{0}\{0}_{1}_rf_model_fold{2}.pkl".format(chembl_id,prot_ft,fold),"wb") as file:
    #     pickle.dump(clf,file)
            
    #test model
    test_pred = clf.predict(X_test)

    accuracy_test = accuracy_score(y_test,test_pred)
    precision_test = precision_score(y_test,test_pred)
    recall_test = recall_score(y_test,test_pred)
    f1_test = f1_score(y_test,test_pred)  
    mcc_test = matthews_corrcoef(y_test,test_pred)

    # #save predictions
    # pred_file = test_ft.iloc[:,:2]
    # pred_file["predicted_bioactivity"] = test_pred
    # pred_file.to_csv(r"..\predictions\small_scale\{0}_{1}_rf_predictions_fold{2}.tsv".format(chembl_id,prot_ft,fold),sep="\t",index=None)

    end = datetime.now()
    print("elapsed time: {}".format(end-start))    
    return [accuracy_test,precision_test,recall_test,f1_test,mcc_test]
    

print(datetime.now())

def CompCentric_model_generation(folder,chembl_id):
    protein_features = ["aac","aac_pssm","aadp_pssm","aatp_pssm","ab_pssm","apaac","cksaagp",\
                        "cksaap","ctdc","ctdd","ctdt","ctriad","dde","dpc","dpc_pssm",\
                            "dp_pssm","d_fpssm","edp_pssm","eedp_pssm","gaac","gdpc","geary","gtpc",\
                                "ksctriad","k-sep_pssm","medp_pssm","moran","nmbroto","paac",\
                                    "pfam","pse_pssm","pssm_ac","pssm_cc","pssm_composition",\
                                        "qso","random200","rpm_pssm","rpssm","spmap","taap",\
                                            "tpc","tpc_pssm","tri-gram_pssm"]           

    test_results_mean = []
    for pft in protein_features: 
        test_results = []
        for i in range(1,6):
            train = pd.read_csv(r"..\datasets\small_scale\{0}\{1}_train_fold{2}.tsv".format(folder,chembl_id,i),sep="\t")
            test = pd.read_csv(r"..\datasets\small_scale\{0}\{1}_test_fold{2}.tsv".format(folder,chembl_id,i),sep="\t")
    
            model = RF_model(train,test,pft,chembl_id,i)       
            test_results.append(model)
        test_results_mean.append([pft]+list(np.mean(test_results, axis=0)))

    df_test_results = pd.DataFrame(columns=["model","accuracy","precision","recall","f1-score","MCC"],data=test_results_mean)            
    df_test_results.to_csv(r"..\results\small-scale_{0}_rf_test_results.tsv".format(chembl_id),sep="\t",index=None)
        

  
chembl44 = CompCentric_model_generation("ChEMBL44_Genistein", "chembl44")
chembl50 = CompCentric_model_generation("ChEMBL50_Quercetin", "chembl50")
chembl83 = CompCentric_model_generation("ChEMBL83_Tamoxifen", "chembl83")
chembl91 = CompCentric_model_generation("ChEMBL91_Miconazole", "chembl91")
chembl104 = CompCentric_model_generation("ChEMBL104_Clotrimazole", "chembl104")
chembl633 = CompCentric_model_generation("ChEMBL633_Amiodarone", "chembl633")
chembl808 = CompCentric_model_generation("ChEMBL808_Econazole", "chembl808")
chembl116438 = CompCentric_model_generation("ChEMBL116438_Curcumin", "chembl116438")
chembl295698 = CompCentric_model_generation("ChEMBL295698_Levoketoconazole", "chembl295698")


print(datetime.now())

















