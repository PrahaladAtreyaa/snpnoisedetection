def testData(bed_file_path,ignore_regions_file_path,germline_file_path,somatic_file_path,bam_file_path,model_file_path):
    
    from GenUtils import getNegRegions
    from GenUtils import getNegDataFrame
    from GenUtils import getModelStats
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_auc_score,roc_curve
    
    dfPosNonGerm=getPosDataFrame(somatic_file_path,bam_file_path,is_germline=False)
    dfPosGerm=getPosDataFrame(germline_file_path,bam_file_path,True)
    dfPos=pd.concat([dfPosGerm,dfPosNonGerm],ignore_index=True)
    
    neg_regions=getNegRegions(bed_file_path,ignore_regions_file_path,germline_file_path,somatic_file_path)
    dfNeg=getNegDataFrame(neg_regions,bam_file_path)
    
    dfFull=pd.concat([dfNeg,dfPos],ignore_index=True)
    # load the model from disk
    rnd = joblib.load(model_file_path)
    dfFull=preprocessTrain(dfFull)
    
    unique_features=['QNAME','REGION','LABEL']
    
    X=dfFull.drop(unique_features,axis=1)
    y=dfFull.LABEL
    y=y.astype('int')
    
    pred=rnd.predict_proba(X)[:, 1]
    dfFull['PROB']=pred
    
    
    #Testing Data Stats
    getScoreStats(rnd,X,y)
    `
    fpr, tpr, thresholds = roc_curve(y.ravel(),pred, pos_label=1)
    auc = roc_auc_score(y.ravel(),pred)
    plt.plot(fpr,tpr,label="auc="+str(auc))
    plt.legend(loc=4)
    plt.show()
    