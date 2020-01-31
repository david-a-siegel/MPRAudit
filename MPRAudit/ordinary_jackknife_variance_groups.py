def ordinary_jackknife_variance_groups(list_of_values_0,list_of_values_1):
    list_of_values_0 = np.array(list_of_values_0)
    list_of_values_1 = np.array(list_of_values_1)
    N = len(list_of_values_0)+len(list_of_values_1)    

    x_i = []
    for i in range(len(list_of_values_0)):
        list_of_values_i = np.append(list_of_values_0[:i],list_of_values_0[i+1:])
        x_i.append((len(list_of_values_i)*np.mean(list_of_values_i)+len(list_of_values_1)*np.mean(list_of_values_1))/(len(list_of_values_i)+len(list_of_values_1)))
    for i in range(len(list_of_values_1)):
        list_of_values_i = np.append(list_of_values_1[:i],list_of_values_1[i+1:])
        x_i.append((len(list_of_values_i)*np.mean(list_of_values_i)+len(list_of_values_0)*np.mean(list_of_values_0))/(len(list_of_values_i)+len(list_of_values_0)))
    
    x_i = np.array(x_i)
    jack_var_x = float(N-1)/N*sum((x_i-np.mean(x_i))**2)
    return jack_var_x
