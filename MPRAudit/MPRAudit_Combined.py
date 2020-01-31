
#Need to create a function that returns R/(R+D) and is editable by the user.

def ratio_function(array1,array2):
    RNA_counts = np.array(array1)
    DNA_counts = np.array(array2)
    if len(RNA_counts)!=len(DNA_counts):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    return RNA_counts/(RNA_counts+DNA_counts)
    
    
def MPRAudit_Pairs(RNA_counts1, DNA_counts1, RNA_counts2, DNA_counts2, sequence_indicators) #sequence_indicator groups clones into otherwise identical sequences
    if len(RNA_counts1)!=len(DNA_counts1):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    if len(RNA_counts1)!=len(RNA_counts2):
        raise Exception('RNA_counts1 and RNA_counts2 must be array-like of the same length')
    if len(RNA_counts1)!=len(DNA_counts2):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')        
    if len(sequence_indicators)!=len(RNA_counts):
        raise Exception('sequence_indicators must be array-like of the same length as the counts')

    RNA_counts1 = np.array(RNA_counts1)
    DNA_counts1 = np.array(DNA_counts1)
    RNA_counts2 = np.array(RNA_counts2)
    DNA_counts2 = np.array(DNA_counts2)
    sequence_indicators = np.array(sequence_indicators)
    
    
    jackknife_variance_list = []
    RD_list = []
    
    for sequence in np.unique(sequence_indicators):
        RNA_sequence_counts1 = RNA_counts1[sequence_indicators==sequence]
        DNA_sequence_counts1 = DNA_counts1[sequence_indicators==sequence]    
        RNA_sequence_counts2 = RNA_counts2[sequence_indicators==sequence]
        DNA_sequence_counts2 = DNA_counts2[sequence_indicators==sequence]    
        jackknife_variance_list.append(deleteDjackknife_variance_T0(RNA_sequence_counts1, DNA_sequence_counts1)+deleteDjackknife_variance_T0(RNA_sequence_counts2, DNA_sequence_counts2))
        RD_list.append(float(sum(RNA_sequence_counts1))/(sum(RNA_sequence_counts1)+sum(DNA_sequence_counts1))-float(sum(RNA_sequence_counts2))/(sum(RNA_sequence_counts2)+sum(DNA_sequence_counts2)))
        
    tech_variance = np.mean(jackknife_variance_list)
    total_variance = np.var(RD_list)
    delta_tech = ordinary_jackknife_variance(jackknife_variance_list)
    delta_var = ordinary_variance_jackknife_variance(RD_list))
    b2_mean, b2_var = uncertainty(tech_variance,total_variance,delta_tech,delta_var)
    return b2_mean, np.sqrt(b2_var)
    
    
def MPRAudit_Groups(RNA_counts, DNA_counts, sequence_indicators, group_indicators=None) #sequence_indicator groups clones into otherwise identical sequences
    if len(RNA_counts)!=len(DNA_counts):
        raise Exception('RNA_counts and DNA_counts must be array-like of the same length')
    if len(sequence_indicators)!=len(RNA_counts):
        raise Exception('sequence_indicators must be array-like of the same length as the counts')
    if group_indicators is not None:
        if len(group_indicators)!=len(RNA_counts):
            raise Exception('sequence_indicators must be array-like of the same length as the counts')    

    RNA_counts = np.array(RNA_counts)
    DNA_counts = np.array(DNA_counts)
    sequence_indicators = np.array(sequence_indicators)
    group_indicators = np.array(group_indicators)
    #Check lengths?
    
    tech_variance = 0
    total_variance = 0
    for group in np.unique(group_indicators):
        jackknife_variance_list = []
        RD_list = []
        for sequence in np.unique(sequence_indicators):
            RNA_sequence_counts = RNA_counts[(sequence_indicators==sequence)&(group_indicators==group)]
            DNA_sequence_counts = DNA_counts[(sequence_indicators==sequence)&(group_indicators==group)]
            jackknife_variance_list.append(deleteDjackknife_variance_T0(RNA_sequence_counts, DNA_sequence_counts))
            RD_list.append(float(sum(RNA_sequence_counts))/(sum(RNA_sequence_counts)+sum(DNA_sequence_counts)))
        tech_variance += np.mean(jackknife_variance_list)*len(jackknife_variance_list)
        total_variance += np.var(RD_list)*len(RD_list)
    tech_variance = tech_variance/len(RNA_counts)
    total_variance = total_variance/len(RNA_counts)
    delta_tech = ordinary_jackknife_variance_Ngroups(jackknife_variance_list,group_indicators)
    delta_var = ordinary_variance_jackknife_variance_Ngroups(RD_list,group_indicators)
    b2_mean, b2_var = uncertainty(tech_variance,total_variance,delta_tech,delta_var)
    return b2_mean, np.sqrt(b2_var)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    