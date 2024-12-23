def getGCScores(samfile,pos,region):
    
    gc_scores=[]
    for i,read in enumerate(samfile.fetch(region=region)):
        chars=[]
        
        if(read.is_reverse):
            for index in range(1,11):
                if(((pos[i]-index)<=-1)==False):
                    chars.append([el for el in reversed(read.get_forward_sequence())][pos[i]-index])
                else:
                    break
            chars=complement(chars)
        else:
            for index in range(1,11):
                if(((pos[i]-index)<=-1)==False):
                    chars.append(read.get_forward_sequence()[pos[i]-index])
                else:
                    break

        gc_scores.append(getGCScore(chars))
        
    return gc_scores

def extractBaseQual(df_comb_features):
    
    base_quality=[]
    for i in range(df_comb_features.shape[0]):
        index=df_comb_features['B_POS'].iloc[i]
        base_quality.append([el for el in reversed(df_comb_features['BASE_QUALITIES'].iloc[i])][index])
    
    df_comb_features['BASE_QUALITY']=base_quality
    
    return df_comb_features

def getGCScore(chars):
    
    C_count = 0
    G_count = 0
    total = len(chars)
    gc_score=0
    for i in chars:
        
        if i=='C':
            C_count+=1
        elif i=='G':
            G_count+=1
    if(total!=0):
        gc_score = (G_count+C_count)/total
    
    return gc_score

def getContextScore(chars,ch):
    
    score=0
    for c in chars:
        if(c!=ch):
            break
        if(c==ch):
            score+=1
    
    return score

def context(samfile,pos,region):
    
    context_scores=[]
    for i,read in enumerate(samfile.fetch(region=region)):
        
        c1=None
        c2=None
        c3=None
        chars=[]
        ch=read.get_forward_sequence()[pos[i]]
        
        if(((pos[i]-1)<=-1)==False):
            c1=read.get_forward_sequence()[pos[i]-1]
        if(((pos[i]-2)<=-1)==False):
            c2=read.get_forward_sequence()[pos[i]-2]
        if(((pos[i]-3)<=-1)==False):
            c3=read.get_forward_sequence()[pos[i]-3]
        chars=[c1,c2,c3]   
        if(read.is_reverse):
            chars=complement(chars)
        context_scores.append(getContextScore(chars,ch))
        
    return context_scores

def complement(chars):
    
    for index in range(len(chars)):
        if(chars[index]!=None):
            if chars[index].upper()=='C':
                chars[index] = 'G'
            elif chars[index].upper()=='G':
                chars[index] = 'C'
            elif chars[index].upper()=='T':
                chars[index] = 'A'
            elif chars[index].upper()=='A':
                chars[index] = 'T'
            
    return chars

    
def getRefBase(read,region):
    
    return read.get_reference_sequence()[read.get_reference_positions().index(int(region.split('-')[-1])-1)].upper()

def getReadFeatures(samfile,region):

    genome = pysam.Fastafile('hg19.fa')
    ref_base = genome.fetch(region.split(':')[0], int(region.split('-')[-1])-1, int(region.split('-')[-1]))
 
    rows_list = []
    for read in samfile.fetch(region=region):

            dict1 = {}
            dict1['QNAME']=read.query_name
            dict1['MAPQ']=read.mapping_quality
            dict1['CIGAR_STATS']=read.get_cigar_stats()
            dict1['QUERY_LENGTH']=read.infer_query_length()
            dict1['FLAG']=read.flag
            dict1['ALIGN_LENGTH']=read.query_alignment_length
            dict1['AVG_QUERY_ALIGN_QUAL']=sum(np.array(read.query_alignment_qualities)/len(np.array(read.query_alignment_qualities)))
            dict1['IS_SECONDARY']=read.is_secondary
            #dict1['GC_SCORE']=getGCscore(read)
            dict1['REF_BASE']=ref_base.upper()
            dict1['IS_READ1']=read.is_read1
            dict1['FAMILY_SIZE']=read.tags[-1][1]
            dict1['REGION']=region
            dict1['BASE_QUALITIES']=read.query_qualities
            dict1['IS_REVERSE']=read.is_reverse
            #dict1['QUERY_SEQ']=read.get_forward_sequence()
            rows_list.append(dict1)

    df_read_features = pd.DataFrame(rows_list)
    
    return df_read_features,ref_base

def getBaseFeatures(samfile,region):
    
    quality=list()
    seqeunces=list()
    pos=list()
    query_names=list()
    context_score=list()
    ref_base=[]
    gc_score=list()
    
    for pileupcolumn in samfile.pileup(region.split(':')[0], region=region,max_depth = 1000000,stepper='nofilter',min_base_quality=0,min_mapping_quality=0):
    
        if(pileupcolumn.pos==int(region.split('-')[-1])-1):
            print ("\ncoverage at base %s = %s" %
               (pileupcolumn.pos, pileupcolumn.n))
                   
            quality=pileupcolumn.get_query_qualities()
            seqeunces=pileupcolumn.get_query_sequences(mark_matches=True,mark_ends=True,add_indels=True)
            pos=pileupcolumn.get_query_positions()
            context_score=context(samfile,pos,region)
            query_names=pileupcolumn.get_query_names()
            gc_score=getGCScores(samfile,pos,region)
            
            break
    
        
    return (pos,seqeunces,quality,query_names,context_score,gc_score)

def extractFLAG(df_comb_features):
    
    def get_binary(x):
        return '{0:012b}'.format(x)
  
    df_comb_features.FLAG=df_comb_features.FLAG.apply(lambda x:get_binary(x))
    
    df_comb_features['F_SUPP_ALIGN']=df_comb_features.FLAG.apply(lambda x:int(x[0]))
    df_comb_features['F_PROPER_PAIR']=df_comb_features.FLAG.apply(lambda x:int(x[10]))
    df_comb_features['F_NOT_PRIMARY_ALIGN']=df_comb_features.FLAG.apply(lambda x:int(x[3]))
    df_comb_features['F_MATE_UNMAPPED']=df_comb_features.FLAG.apply(lambda x:int(x[8]))
    df_comb_features['F_READ_UNMAPPED']=df_comb_features.FLAG.apply(lambda x:int(x[9]))
    
    df_comb_features.drop(columns=['FLAG'],axis=1,inplace=True)

    return df_comb_features

def extractCIGAR(df_comb_features):
    
    df_comb_features['C_BAM_CMATCH']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][0])
    df_comb_features['C_BAM_CINS']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][1])
    df_comb_features['C_BAM_CDEL']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][2])
    df_comb_features['C_BAM_CREF_SKIP']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][3])
    df_comb_features['C_BAM_CSOFT_CLIP']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][4])
    df_comb_features['C_BAM_CHARD_CLIP']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][5])
    df_comb_features['C_BAM_CPAD']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][6])
    df_comb_features['C_BAM_CEQUAL']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][7])
    df_comb_features['C_BAM_CDIFF']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][8])
    df_comb_features['C_BAM_CBACK']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][9])
    #df_comb_features['C_NM_TAG']=df_comb_features.CIGAR_STATS.apply(lambda x:x[0][10])
    
    df_comb_features.drop(columns=['CIGAR_STATS'],axis=1,inplace=True)
    
    return df_comb_features

def extractSEQ(x):


    forward_match=False
    reverse_match=False
    forward_mismatch=False
    reverse_mismatch=False
    ref_skip=False
    is_insert=False
    is_del=False
    insert_len=0
    del_len=0
    is_start_of_read=False
    is_end_of_read=False
    
    findex=0
    lindex=len(x)

    
    if('-' in x):
        
        last_index=x.index('-')
        
    if('+' in x):
        
        last_index=x.index('+')
        is_insert=True
        
    if(x.startswith('^')):
        
        first_index=2
        is_start_of_read=True
        
    if('*' in x):
        
        is_del=True
    
    if(x.endswith('$')):
        
        is_end_of_read=True
    
 
    if ('a' in x[findex:lindex] or 'g' in x[findex:lindex] or 'c' in x[findex:lindex]  or 't' in x[findex:lindex]):  
        reverse_mismatch=True
    if ('A' in x[findex:lindex] or 'G' in x[findex:lindex] or 'C' in x[findex:lindex] or 'T' in x[findex:lindex]): 
        forward_mismatch=True
        
    return (forward_match,reverse_match,forward_mismatch,reverse_mismatch,ref_skip,is_insert,is_del,insert_len,del_len,is_start_of_read,is_end_of_read)

def getBaseRelPos(row):
    
    query_length=row['QUERY_LENGTH']
    base_pos=row['B_POS']
    
    return base_pos/query_length

def labelSubstitutionDeletion(reference,supporting_variant,seq):
    
    index=None
    label=None
    
    if(len(seq)==1 or '-' in seq[1] or '+' in seq[1]):
        index=0
    elif('-' in seq):
        index=seq.index('-')-1
    elif('+' in seq):
        index=seq.index('+')-1
    elif(seq.startswith('^')):
        index=2
    elif(seq.endswith('$')):
        index=-2
    else:
        raise ValueError('Wrong sequence:'+seq)
        
    if(seq[index].upper()==reference.upper()):
        label="M"
    elif(seq[index].upper() in supporting_variant):
        label="SVR"
    else:
        label="NSVR"
        
    return label

def labelInsertion(reference,supporting_variant,seq):
    
    label=None
    ref_index=None
    ins_index=None
    
    if('+' not in seq):
        label="M"
    elif('+' in seq):
        ref_index=seq.index('+')-1
        ins_index=seq.index('+')+2
        if(seq[ins_index].upper() in supporting_variant):
            label="SVR"
        else:
            label="NSVR"
    return label    

def negativeLabels(reference,seq):
    
    index=None
    label=None
    #print(seq)
    if(len(seq)==1 and seq[0]=="^"):
        return "ERR"
    if(len(seq)==1 and seq[0]=="$"):
        return "ERR"
    
    if('+' in seq or '*' in seq):
        label="NSVR"
        return label
    
    if(len(seq)==1 or '-' in seq[1]):
        index=0
    elif('-' in seq):
        index=seq.index('-')-1
    elif(seq.startswith('^')):
        index=2
    elif(seq.endswith('$')):
        index=-2
    else:
        raise ValueError('Wrong sequence:'+seq)
     
    if(seq[index].upper()==reference.upper()):
        label="M"
    else:
        label="NSVR"
        
    return label

def getDataframe(file_name,region,train_mode=True,variant_type=None,supporting_variant=None,is_control_positive=True):
    

    samfile = pysam.AlignmentFile(file_name, "rb")
    df_read_features,reference=getReadFeatures(samfile,region)
    df_comb_features=df_read_features
    
    (pos,seqeunces,quality,query_names,context_score,gc_score)=getBaseFeatures(samfile,region)
    df_comb_features['B_POS']=pos
    df_comb_features['B_SEQ']=[c for c in seqeunces if c!="^"]
    df_comb_features['B_BASE_QUAL']=quality
    df_comb_features['B_BQNAME']=query_names
    df_comb_features['B_CONTEXT_SCORE']=context_score
    df_comb_features['B_GC_SCORE']=gc_score
    
    df_comb_features.drop(columns=['B_BQNAME'],axis=1,inplace=True)
    
    df_comb_features=extractBaseQual(df_comb_features)
    df_comb_features=extractFLAG(df_comb_features)
    
    
    df_comb_features=extractCIGAR(df_comb_features)
    

    df_comb_features['S_FWD_MISMATCH']=df_comb_features.B_SEQ.apply(lambda x:extractSEQ(x)[2])
    df_comb_features['S_REV_MISMATCH']=df_comb_features.B_SEQ.apply(lambda x:extractSEQ(x)[3])
    df_comb_features['S_IS_INS']=df_comb_features.B_SEQ.apply(lambda x:extractSEQ(x)[5])
    df_comb_features['S_IS_DEL']=df_comb_features.B_SEQ.apply(lambda x:extractSEQ(x)[6])
    df_comb_features['S_IS_START']=df_comb_features.B_SEQ.apply(lambda x:extractSEQ(x)[9])
    df_comb_features['S_IS_END']=df_comb_features.B_SEQ.apply(lambda x:extractSEQ(x)[10])

    
    df_comb_features['B_REL_POS']=df_comb_features.apply(lambda row: getBaseRelPos(row),axis=1)
    
    if(train_mode==True):
        if(is_control_positive):
            if(variant_type=="Substitution" or variant_type=="Deletion"):
                df_comb_features['LABEL']=df_comb_features.B_SEQ.apply(lambda x:labelSubstitutionDeletion(reference,supporting_variant,x))
            else:

                df_comb_features['LABEL']=df_comb_features.B_SEQ.apply(lambda x:labelInsertion(reference,supporting_variant,x))
        else:

            df_comb_features['LABEL']=df_comb_features.B_SEQ.apply(lambda x:negativeLabels(reference,x))
    
 
    #df_comb_features.drop(columns=['QUERY_LENGTH','ALIGN_LENGTH','B_POS','B_SEQ'],axis=1,inplace=True)
    
    df_comb_features.drop(columns=['BASE_QUALITIES'],axis=1,inplace=True)
    
    return df_comb_features


def bed_to_region(bed_file):
    
    bed_regions=[]
    for row in range(len(bed_file)):

        chromosome=bed_file.iloc[row]['Chromosome']
        start=bed_file.iloc[row]['Start']
        end=bed_file.iloc[row]['End']

        for reg in range(int(start),int(end)+1):
            bed_region=chromosome+":"+str(reg)+"-"+str(reg)
            bed_regions.append((bed_region))
            
    return bed_regions

def df_to_region(df):
    
    region_list=[]
    
    for row in range(len(df)):
    
        chromosome=df.iloc[row]['Chromosome']
        start=df.iloc[row]['Start']
        end=df.iloc[row]['End']
        region_list.append(chromosome+":"+str(start)+"-"+str(end))
        
    return region_list

def val_to_region(file_name,is_germline):
    
    from Fileutils import multi_to_single_base_deletion

    df=multi_to_single_base_deletion(file_name,is_germline)
    df=df_to_region(df)
    
    return df

#Util function to extract a list of regions where no variation is expected from all aligned regions using
#the formula all_aligned_regions-(union(germline_regions,somatic_regions))-ignore_regions

def getNegRegions(bed_file_path,ignore_regions_file_path,germline_file_path,somatic_file_path):
    
    bed_file=pd.read_csv(bed_file_path,sep='\t',header=None,names=['Chromosome','Start','End','NA'])
    bed_regions=bed_to_region(bed_file)
    unique_regions_file=pd.read_csv(ignore_regions_file_path,sep='\t')
    ignore_regions=df_to_region(unique_regions_file)
    pos_regions_nongerm=val_to_region(somatic_file_path,False)
    pos_regions_germ=val_to_region(germline_file_path,True)
    pos_regions=list(set(pos_regions_germ).union(set(pos_regions_nongerm)))
    pos_regions= sorted(pos_regions, key = lambda x: x.split(':')[0])
    neg_regions=set(bed_regions)-set(ignore_regions)-set(pos_regions)
    neg_regions= sorted(neg_regions, key = lambda x: x.split(':')[0])
    
    return neg_regions


#To get Data from negative control regions where no variations are expected
def getNegDataFrame(regions,bam_name):
    
    #neg_file=pd.read_csv(panel_name,sep='\t',skiprows=2)
    
    dataframes_n=[]
    counter=0
    for index,region in enumerate(regions):
        dfFull=getDataframe(bam_name,region=region, train_mode=True,is_control_positive=False)
        dfPositive=dfFull[dfFull.LABEL=="M"]
        dfNegative=dfFull[dfFull.LABEL!="M"]
        
        
        df=pd.concat([dfPositive.sample(n=int(dfPositive.shape[0]/5), random_state=0),dfNegative],ignore_index=True)
        counter=counter+df.shape[0]
        print("=================")
        print(dfNegative.shape[0])
        print(index)
        print(region)
        print(df.shape[0])
        print(counter)
        print("================")
        dataframes_n.append(df)
    df=pd.concat(dataframes_n,ignore_index=True)
    
    return df

def preprocessTrain(df):
    
    import pandas as pdd

    df=df.drop(['F_SUPP_ALIGN','F_PROPER_PAIR','F_MATE_UNMAPPED','F_READ_UNMAPPED',
            'C_BAM_CREF_SKIP','C_BAM_CHARD_CLIP','C_BAM_CPAD',
            'C_BAM_CEQUAL','C_BAM_CDIFF','C_BAM_CBACK','F_NOT_PRIMARY_ALIGN'],axis=1)

    df=pdd.get_dummies(df,columns=['REF_BASE'])
    bases=['A','G','C','T']
    
    for base in bases:
        if("REF_BASE_"+base not in df.columns):
            df["REF_BASE_"+base]=0
    
    imp_features=['QNAME','REGION','B_POS','ALIGN_LENGTH', 'AVG_QUERY_ALIGN_QUAL', 'FAMILY_SIZE', 'IS_READ1',
       'IS_REVERSE', 'MAPQ', 'QUERY_LENGTH', 'B_CONTEXT_SCORE','C_BAM_CSOFT_CLIP',
       'B_GC_SCORE', 'BASE_QUALITY', 'C_BAM_CMATCH', 
       'C_BAM_CDEL', 'S_FWD_MISMATCH', 'S_REV_MISMATCH',
       'S_IS_DEL', 'B_REL_POS', 'REF_BASE_A',
       'REF_BASE_C', 'REF_BASE_G', 'REF_BASE_T','LABEL']
    df['LABEL'][df.LABEL=="M"]=1
    df['LABEL'][df.LABEL=="SVR"]=1
    df['LABEL'][df.LABEL=="NSVR"]=0
    
    return df[imp_features]

def getModelStats(estimator,X,y):
    
    from sklearn.metrics import confusion_matrix,auc,roc_curve,accuracy_score
    
    print("------------------------------------")
    tn, fp, fn, tp = confusion_matrix(y, estimator.predict(X)).ravel()
    print("fpr="+str(fp/(tn+fp)))
    print("acc="+str(accuracy_score(y,estimator.predict(X))))
    fpr, tpr, _ = roc_curve(y, estimator.predict_proba(X)[:, 1], pos_label=1)
    print("auc="+str(auc(fpr, tpr)))
    print("tpr="+str(tp/(fn+tp)))
    print("------------------------------------")
    print("Confusion Matrix")
    print("tp="+str(tp))
    print("tn="+str(tn))
    print("fp="+str(fp))
    print("fn="+str(fn))     
    print("------------------------------------")

def getNegDataFrame(regions,bam_name):
    
    #neg_file=pd.read_csv(panel_name,sep='\t',skiprows=2)
    
    dataframes_n=[]
    counter=0
    for index,region in enumerate(regions):
        dfFull=getDataframe(bam_name,region=region, train_mode=True,is_control_positive=False)
        dfPositive=dfFull[dfFull.LABEL=="M"]
        dfNegative=dfFull[dfFull.LABEL!="M"]
        
        
        df=pd.concat([dfPositive,dfNegative],ignore_index=True)
        counter=counter+df.shape[0]
        dataframes_n.append(df)
    df=pd.concat(dataframes_n,ignore_index=True)
    
    return df