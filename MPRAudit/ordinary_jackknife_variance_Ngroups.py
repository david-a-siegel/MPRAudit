def ordinary_jackknife_variance_Ngroups(list_of_values,group_indicators):
    list_of_values = np.array(list_of_values)
    group_indicators = np.array(group_indicators)
    N = len(list_of_values)
    
    x_i = []
    for group in np.unique(group_indicators):
        group_values = list_of_values[group_indicators==group]
        notgroup_values = list_of_values[group_indicators!=group]
        for i in range(len(group_values)):
            group_values_i = np.append(group_values[:i],group_values[i+1:])
            x_i.append(np.mean(np.append(group_values_i,notgroup_values)))
    
    x_i = np.array(x_i)
    jack_var_x = float(N-1)/N*sum((x_i-np.mean(x_i))**2)
    return jack_var_x


