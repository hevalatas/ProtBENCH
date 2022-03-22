# -*- coding: utf-8 -*-
import pandas as pd
# import pickle

from sklearn.ensemble import RandomForestRegressor
from score_metrics import rmse,spearman,f1,mcc
from datetime import datetime



def PCM_model(tr,ts,prot_ft,comp_ft="ecfp4",pKd=7.0):
    start = datetime.now()
    
    df_prot_ft = pd.read_csv(r"..\datasets\medium_scale\feature_vectors\{}.tsv".format(prot_ft),sep="\t")
    df_comp_ft = pd.read_csv(r"..\datasets\medium_scale\feature_vectors\{}.tsv".format(comp_ft),sep="\t")

    train_ft = tr.merge(df_prot_ft,on="target_id").merge(df_comp_ft,on="compound_id")
    X_train = train_ft.iloc[:,3:]
    y_train = train_ft["-log[M]"]
    
    test_ft = ts.merge(df_prot_ft,on="target_id").merge(df_comp_ft,on="compound_id")
    X_test = test_ft.iloc[:,3:]
    y_test = test_ft["-log[M]"].values
   
    #estimator
    reg = RandomForestRegressor(n_estimators=100,max_features=0.33,random_state=42) 
    reg.fit(X_train, y_train)
    test_pred = reg.predict(X_test)

    model_name = prot_ft if comp_ft=="ecfp4" else prot_ft+"_"+comp_ft   
    
    #save model
    # with open(r"..\models\medium_scale\{}_rf.pkl".format(model_name),"wb") as file:
    #     pickle.dump(reg,file)

    rmse_test = rmse(y_test,test_pred)
    spearman_test = spearman(y_test,test_pred)
    f1_test = f1(y_test,test_pred,pKd)    
    mcc_test = mcc(y_test,test_pred,pKd)

    # #save predictions
    # pred_file = test_ft.iloc[:,:3]
    # pred_file["predicted_value"] = test_pred
    # pred_file.to_csv(r"..\predictions\medium_scale\{}_predictions.tsv".format(model_name),sep="\t",index=None)

    end = datetime.now()
    print("elapsed time: {}".format(end-start))    
    return [model_name,rmse_test,spearman_test,f1_test,mcc_test]


train = pd.read_csv(r"..\datasets\medium_scale\mDavis_train.tsv",sep="\t")
test = pd.read_csv(r"..\datasets\medium_scale\mDavis_test.tsv",sep="\t")

protein_features = ["apaac","ctdd","ctriad","dde","geary","k-sep_pssm","pfam","qso","random200","spmap","taap",\
                  "protvec","seqvec","transformer-avg","transformer-pool","unirep1900","unirep5700"]           

test_results = []
for pft in protein_features:
    model = PCM_model(train,test,pft)    
    test_results.append(model)
            
df_test_results = pd.DataFrame(columns=["model","RMSE","Spearman","F1-score","MCC"],data=test_results)
random_model = PCM_model(train,test,"random200","random-ecfp4")
df_test_results.loc[len(df_test_results)] = random_model
df_test_results.to_csv(r"..\results\medium_scale\mDavis_rf-models_test-results.tsv",sep="\t",index=None)
  
 










            




