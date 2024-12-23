
#To extract information of reference base and expected variant from control files of both germaline and somatic varaints to label data

def getPosDataFrame(panel_name,bam_name,is_germline):

    dataframes_p=[]
    
    pos_file=multi_to_single_base_deletion(panel_name,is_germline)

    for row in range(pos_file.shape[0]):
        chromosome=pos_file.iloc[row].Chromosome
        ref_pos=pos_file.iloc[row].End
        region=chromosome+":"+str(ref_pos)+"-"+str(ref_pos)
        variant_type=pos_file.iloc[row]['Variant Type']

        #reference=None
        supporting_variant=None

        if(variant_type=="Substitution"):
            #reference=pos_file.iloc[row].Reference
            if(is_germline==False):
                supporting_variant=pos_file.iloc[row]['Variant Allele'].split('/')
            else:
                supporting_variant=pos_file.iloc[row]['SNP Call (HG002)'].split('/')

        elif(variant_type=='Deletion'):
            #reference=pos_file.iloc[row].Reference[1]
            supporting_variant="*".split('/')

        else:
            #reference=pos_file.iloc[row].Reference
            supporting_variant=pos_file.iloc[row]['Variant Allele'][1].split('/')
        
        dataframes_p.append(getDataframe(bam_name,region=region,train_mode=True,variant_type=variant_type,supporting_variant=supporting_variant))

    df=pd.concat(dataframes_p,ignore_index=True)
    

#Util function to convert multi base deletion to single base deletions

def multi_to_single_base_deletion(val_file_name,is_germline):
    
    file=pd.read_csv(val_file_name,sep='\t')
    row_list=[]
    rows_to_drop=[]
    for row in range(file.shape[0]):
        if(file.iloc[row]['Variant Type']=="Deletion" and len(file.iloc[row]['Reference'])==2):
            file.iloc[row]['Start']=file.iloc[row]['Start']+1
        if(file.iloc[row]['Variant Type']=="Deletion" and len(file.iloc[row]['Reference'])>2):
            
            
            #print(row)
            context=file.iloc[row]['Reference'][0]
            for index,char in enumerate(file.iloc[row]['Reference'][1:]):
                
                dict1={}
                dict1['Chromosome']=file.iloc[row]['Chromosome']
                dict1['Start']=file.iloc[row]['Start']+index
                dict1['End']=file.iloc[row]['Start']+index+1
                dict1['Reference']=context+char
                #print(file.iloc[row]['Variant Allele'])
                dict1['Variant Allele']=file.iloc[row]['Variant Allele']
                dict1['Variant Type']=file.iloc[row]['Variant Type']
                if(is_germline):
                    dict1['SNP Call (HG002)']=file.iloc[row]['SNP Call (HG002)']
                    dict1['Score (HG002)']=file.iloc[row]['Score (HG002)']
                    dict1['Zygosity (HG002)']=file.iloc[row]['Zygosity (HG002)']
                row_list.append(dict1)
            rows_to_drop.append(row)
    file=file.drop(rows_to_drop)
    df = pd.DataFrame(row_list)
    df_out=pd.concat([file,df],ignore_index=True)
    
    for row in range(df_out.shape[0]):
        if(df_out.iloc[row]['Variant Type']=="Deletion" and len(df_out.iloc[row]['Reference'])==2):
            df_out['Start'].iloc[row]=df_out['Start'].iloc[row]+1
    
    return df_out

#