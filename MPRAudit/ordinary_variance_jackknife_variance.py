import numpy as np

def ordinary_variance_jackknife_variance(list_of_values):
    list_of_values = np.array(list_of_values)
    
    x_i = []
    for i in range(len(list_of_values)):
        list_of_values_i = np.append(list_of_values[:i],list_of_values[i+1:])
        x_i.append(np.var(list_of_values_i))
    x_i = np.array(x_i)
    jack_var_x = float(len(list_of_values)-1)/len(list_of_values)*sum((x_i-np.mean(x_i))**2)
    return jack_var_x