def deleteDjackknife_variance_T0(RNA_counts, DNA_counts):
    if len(RNA_counts)!=len(DNA_counts):
        print("Error: RNA counts has different length from DNA counts")
    RNA_counts = np.array(RNA_counts)
    DNA_counts = np.array(DNA_counts)
    
    num_trials = 100
    
    exp_pow = 3./5
    
    d = int(len(RNA_counts)**(exp_pow))
    d_kept = len(RNA_counts) - d #This is n-d
    
    x_i = []
    for i in range(num_trials):
        kept_clones = random.sample(range(len(RNA_counts)),d_kept)
        RNA_counts_i = RNA_counts[kept_clones]
        DNA_counts_i = DNA_counts[kept_clones]
        x_i.append(float(sum(RNA_counts_i))/(sum(RNA_counts_i)+sum(DNA_counts_i)))
    x_i = np.array(x_i)
    jack_var_x = float(d_kept)/d/num_trials*sum((x_i-np.mean(x_i))**2)
    return jack_var_x

