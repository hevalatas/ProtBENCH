# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from datetime import datetime
# import pickle

from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import MinMaxScaler
from score_metrics import rmse,spearman,mcc, multiclass_mcc


print(datetime.now())

def PCM_model(tr,ts,tr_pchembl_median,prot_ft,comp_ft="ecfp4"):
    start = datetime.now()

    df_comp_ft = pd.read_csv(r"..\datasets\large_scale\feature_vectors\{}.tsv".format(comp_ft),sep="\t")
    tr_comp_ft = df_comp_ft.loc[df_comp_ft["compound_id"].isin(tr["compound_id"])]
    ts_comp_ft = df_comp_ft.loc[df_comp_ft["compound_id"].isin(ts["compound_id"])]

    if prot_ft == "null":
        train_ft = tr.merge(df_comp_ft,on="compound_id")
        test_ft = ts.merge(df_comp_ft,on="compound_id")
    else:        
        df_prot_ft = pd.read_csv(r"..\datasets\large_scale\feature_vectors\{}.tsv".format(prot_ft),sep="\t")
        tr_prot_ft = df_prot_ft.loc[df_prot_ft["target_id"].isin(tr["target_id"])].copy()
        ts_prot_ft = df_prot_ft.loc[df_prot_ft["target_id"].isin(ts["target_id"])].copy()

        #scale protein feature vectors based on train data
        scaler = MinMaxScaler()   
        tr_prot_ft.iloc[:,1:] = scaler.fit_transform(tr_prot_ft.iloc[:,1:]) 
        ts_prot_ft.iloc[:,1:] = scaler.transform(ts_prot_ft.iloc[:,1:])  
 
        train_ft = tr.merge(tr_prot_ft,on="target_id").merge(tr_comp_ft,on="compound_id")
        test_ft = ts.merge(ts_prot_ft,on="target_id").merge(ts_comp_ft,on="compound_id")

    X_train = train_ft.iloc[:,3:]
    y_train = train_ft["pchembl_value"]

    X_test = test_ft.iloc[:,3:]
    y_test = test_ft["pchembl_value"].values
  
    #train model
    reg = RandomForestRegressor(n_estimators=100,max_features=0.33,random_state=42) 
    reg.fit(X_train, y_train)

    model_name = prot_ft if (comp_ft=="ecfp4" and prot_ft != "null") else "only-"+comp_ft if prot_ft=="null" else prot_ft+"_"+comp_ft 

    #save model
    # with open(r"..\models\large_scale\{0}\{1}\{1}_{2}_rf.pkl".format(split,prot_family,model_name),"wb") as file:
    #     pickle.dump(reg,file)

    #test model
    test_pred = reg.predict(X_test)    
    
    med_cor_test_pred = test_pred + tr_pchembl_median - np.median(test_pred)
               
    rmse_test = rmse(y_test,test_pred)
    med_cor_rmse_test = rmse(y_test,med_cor_test_pred)
    spearman_test = spearman(y_test,test_pred)      
    mcc_test = mcc(y_test,test_pred,tr_pchembl_median) 
    med_cor_mcc_test = mcc(y_test,med_cor_test_pred,tr_pchembl_median) 
    multiclass_mcc_test = multiclass_mcc(y_test,test_pred)
   
    # #save predictions
    # pred_file = test_ft.iloc[:,:3]
    # pred_file["predicted_value"] = test_pred
    # pred_file["med_cor_predicted_value"] = med_cor_test_pred
    # pred_file.to_csv(r"..\predictions\large_scale\{0}\{1}\{1}_{2}_predictions.tsv".format(split,prot_family,model_name),sep="\t",index=None)

    end = datetime.now()
    print("elapsed time: {}".format(end-start))
    
    return [model_name,rmse_test,med_cor_rmse_test,spearman_test,\
             mcc_test,med_cor_mcc_test,multiclass_mcc_test]
   

def ProtFamBased_model_generation(split,prot_family):   
    train = pd.read_csv(r"..\datasets\large_scale\{0}\{1}_train.tsv".format(split,prot_family),sep="\t")
    test = pd.read_csv(r"..\datasets\large_scale\{0}\{1}_test.tsv".format(split,prot_family),sep="\t")

    protein_features = ["apaac","ctdd","ctriad","dde","geary","k-sep_pssm","pfam","qso","random200","spmap","taap",\
                  "protvec","seqvec","transformer-avg","transformer-pool","unirep1900","unirep5700"]           

    test_results = []
    for pft in protein_features:
        model = PCM_model(train,test,np.median(train["pchembl_value"]),pft)       
        test_results.append(model)
                    
    df_test_results = pd.DataFrame(columns=["model","rmse","med_cor_rmse","spearman","mcc","med_cor_mcc","multiclass_mcc"],data=test_results)
    
    #baseline models
    random_model = PCM_model(train,test,np.median(train["pchembl_value"]),"random200","random-ecfp4")
    only_ecfp4_model = PCM_model(train,test,np.median(train["pchembl_value"]),"null")
    only_random_ecfp4_model = PCM_model(train,test,np.median(train["pchembl_value"]),"null","random-ecfp4")

    df_baseline_test_results = pd.DataFrame(columns=["model","rmse","med_cor_rmse","spearman","mcc","med_cor_mcc","multiclass_mcc"],\
                                            data=[random_model,only_ecfp4_model,only_random_ecfp4_model])

    df_all_test_results = pd.concat([df_test_results,df_baseline_test_results])
    df_all_test_results.to_csv(r"..\results\large-scale_{0}_{1}_test_results.tsv".format(split,prot_family),sep="\t",index=None)

   

split = "fully_dissimilar_split"
# split = "dissimilar_compound_split"
# split = "random_split"
  
epi_gen = ProtFamBased_model_generation(split,"epigenetic-regulators")   
#hydro = ProtFamBased_model_generation(split,"hydrolases") 
#ion_ch = ProtFamBased_model_generation(split,"ion-channels")   
#memb_rec = ProtFamBased_model_generation(split,"membrane-receptors") 
#oth_enz = ProtFamBased_model_generation(split,"other-enzymes") 
#ox_red = ProtFamBased_model_generation(split,"oxidoreductases") 
#protease = ProtFamBased_model_generation(split,"proteases") 
#trans_fac = ProtFamBased_model_generation(split,"transcription-factors")    
# transferase = ProtFamBased_model_generation(split,"transferases") 
#transporter = ProtFamBased_model_generation(split,"transporters") 



print(datetime.now())







 #-------------------





