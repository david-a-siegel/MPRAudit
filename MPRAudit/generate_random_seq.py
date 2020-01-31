def generate_random_seq(length_of_seq,has_TATACAG):
    seq_vals = np.random.randint(0,4,size=length_of_seq)
    seq_string = ''
    for i in range(length_of_seq):
        if seq_vals[i]==0:
            seq_string+='A'
        elif seq_vals[i]==1:
            seq_string+='C'
        elif seq_vals[i]==2:
            seq_string+='T'
        elif seq_vals[i]==3:
            seq_string+='G'
        else:
            print('error in string generation')
    if has_TATACAG:
        start_position = np.random.randint(0,length_of_seq-len('TATACAG'))
        seq_string = seq_string[:start_position]+'TATACAG'+seq_string[start_position+len('TATACAG'):]
    return seq_string