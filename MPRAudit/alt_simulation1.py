import numpy as np
import scipy.stats as st
import random
import re
import scipy
from scipy.stats import nbinom, poisson
import sys





def main():

    if len(sys.argv) >= 2:
        namenumber = int(sys.argv[1])
    else:
        namenumber = 1

    pretime = time.clock()

    namestring = 'simulation2_'+str(namenumber)+'.txt'
    text_file = open(namestring,"w")
    text_file.write('Num_clones')
    text_file.write(',')
    text_file.write('Fig2_mean_jackknife') #Q2
    text_file.write(',')
    text_file.write('Fig2_var_RD')
    text_file.write(',')
    text_file.write('Fig2_var_jackknife')
    text_file.write(',')
    text_file.write('Fig2_var_var_RD')
    text_file.write(',')
    text_file.write('Fig2_b2_mean')
    text_file.write(',')
    text_file.write('Fig2_b2_var')
    text_file.write(',')
    text_file.write('Fig3EF_a1') #Q5
    text_file.write(',')
    text_file.write('Fig3EF_b1')
    text_file.write(',')
    text_file.write('Fig3EF_da1')
    text_file.write(',')
    text_file.write('Fig3EF_db1')
    text_file.write(',')
    text_file.write('Fig3EF_h21')
    text_file.write(',')
    text_file.write('Fig3EF_dh21')
    text_file.write(',')
    text_file.write('Fig3EF_a2')
    text_file.write(',')
    text_file.write('Fig3EF_b2')
    text_file.write(',')
    text_file.write('Fig3EF_da2')
    text_file.write(',')
    text_file.write('Fig3EF_db2')
    text_file.write(',')
    text_file.write('Fig3EF_h22')
    text_file.write(',')
    text_file.write('Fig3EF_dh22')
    text_file.write(',')
    text_file.write('Fig3BC_a') #Q4
    text_file.write(',')
    text_file.write('Fig3BCb')
    text_file.write(',')
    text_file.write('Fig3BC_da')
    text_file.write(',')
    text_file.write('Fig3BC_db')
    text_file.write(',')
    text_file.write('Fig3BC_h2')
    text_file.write(',')
    text_file.write('Fig3BC_dh2')
    text_file.write(',')
    text_file.write('Fig4BC_a') #Q6 (Q3 was removed)
    text_file.write(',')
    text_file.write('Fig4BC_b')
    text_file.write(',')
    text_file.write('Fig4BC_da')
    text_file.write(',')
    text_file.write('Fig4BC_db')
    text_file.write(',')
    text_file.write('Fig4BC_h2')
    text_file.write(',')
    text_file.write('Fig4BC_dh2')
    text_file.write('\n')

    for num_clones in range(2,21):
        #Simulate Null Oligos:

        num_seqs = 10000
        seqlength = 15
        c_const = -0.5
        b_const = -0.15
        offset_const = 2

        for cluster_repeats in range(2):

            text_file.write(repr(num_clones))
            text_file.write(',')

            num_Ts = []
            num_TATACAGs = []
            A_list = []

            for _ in range(num_seqs/2):
                #Random sequence followed by TATACAG sequence.
                #First, random sequence:
                randseq = generate_random_seq(seqlength,False)
                num_TATACAGs.append(len(re.findall('(?=TATACAG)',randseq)))
                num_Ts.append(len(re.findall('(?=T)',randseq)))
                A_list.append(offset_const+num_TATACAGs[-1]*c_const-num_Ts[-1]*b_const)
                #Now TATACAG sequence:
                TATACAG_start_loc = np.random.randint(0,9) #0-8 inclusive
                randseq = randseq[:TATACAG_start_loc]+'TATACAG'+randseq[TATACAG_start_loc+7:]
                num_TATACAGs.append(len(re.findall('(?=TATACAG)',randseq)))
                num_Ts.append(len(re.findall('(?=T)',randseq)))
                A_list.append(offset_const+num_TATACAGs[-1]*c_const-num_Ts[-1]*b_const)
    
            #Simulate Counts:
            jackknife_variance_list = []
            RD_list = []
            for A_value in A_list:
                RNA_sampled, DNA_sampled, sum_RNA, sum_DNA = return_counts(num_clones=num_clones,dna_per_clone=20,A=A_value,s=0.0,sigma_mu_ratio=30., K=10000000)  #A is the ratio, K is counts
                RNA_sampled = np.array(RNA_sampled,dtype=float)/sum_RNA
                DNA_sampled = np.array(DNA_sampled,dtype=float)/sum_DNA

                jackknife_variance_list.append(deleteDjackknife_variance_T0(RNA_sampled, DNA_sampled))
                RD_list.append(float(sum(RNA_sampled))/(sum(RNA_sampled)+sum(DNA_sampled)))
    
            num_TATACAGs = np.array(num_TATACAGs)
            num_Ts = np.array(num_Ts)
            A_list = np.array(A_list)
            jackknife_variance_list = np.array(jackknife_variance_list)
            RD_list = np.array(RD_list)


            #Fig 2. How much of the total variance is due to sequence?

            b2_mean, b2_var = uncertainty(np.mean(jackknife_variance_list),np.var(RD_list),ordinary_jackknife_variance(jackknife_variance_list),ordinary_variance_jackknife_variance(RD_list))

            text_file.write(repr(np.mean(jackknife_variance_list)))
            text_file.write(',')
            text_file.write(repr(np.var(RD_list)))
            text_file.write(',')
            text_file.write(repr(ordinary_jackknife_variance(jackknife_variance_list)))
            text_file.write(',')
            text_file.write(repr(ordinary_variance_jackknife_variance(RD_list)))
            text_file.write(',')
            text_file.write(repr(b2_mean))
            text_file.write(',')
            text_file.write(repr(b2_var))
            text_file.write(',')



            #Fig 3EF Do sequences with TATACAG have more b^2 than randomly generated sequences?


            jackknife_variance_list_TATACAG = jackknife_variance_list[np.arange(0,num_seqs,2)+1]
            jackknife_variance_list_noTATACAG = jackknife_variance_list[np.arange(0,num_seqs,2)]
            RD_list_TATACAG = RD_list[np.arange(0,num_seqs,2)+1]
            RD_list_noTATACAG = RD_list[np.arange(0,num_seqs,2)]

            a = np.mean(jackknife_variance_list_TATACAG)
            b = np.var(RD_list_TATACAG)
            da = ordinary_jackknife_variance(jackknife_variance_list_TATACAG)
            db = ordinary_variance_jackknife_variance(RD_list_TATACAG)

            uncert1,uncert2 = uncertainty(a,b,da,db)

            text_file.write(repr(a))
            text_file.write(',')
            text_file.write(repr(b))
            text_file.write(',')
            text_file.write(repr(da))
            text_file.write(',')
            text_file.write(repr(db))
            text_file.write(',')
            text_file.write(repr(uncert1))
            text_file.write(',')
            text_file.write(repr(uncert2))
            text_file.write(',')
    

            a = np.mean(jackknife_variance_list_noTATACAG)
            b = np.var(RD_list_noTATACAG)
            da = ordinary_jackknife_variance(jackknife_variance_list_noTATACAG)
            db = ordinary_variance_jackknife_variance(RD_list_noTATACAG)


            uncert1,uncert2 = uncertainty(a,b,da,db)


            text_file.write(repr(a))
            text_file.write(',')
            text_file.write(repr(b))
            text_file.write(',')
            text_file.write(repr(da))
            text_file.write(',')
            text_file.write(repr(db))
            text_file.write(',')
            text_file.write(repr(uncert1))
            text_file.write(',')
            text_file.write(repr(uncert2))
            text_file.write(',')
    


            #Fig 3BC. If we only know about TATACAG, can we infer some residual b2?

            jackknife_variance_list_0 = jackknife_variance_list[num_TATACAGs==0]
            jackknife_variance_list_1 = jackknife_variance_list[num_TATACAGs==1]
            RD_list_0 = RD_list[num_TATACAGs==0]
            RD_list_1 = RD_list[num_TATACAGs==1]



            a = (len(jackknife_variance_list_0)*np.mean(jackknife_variance_list_0)+len(jackknife_variance_list_1)*np.mean(jackknife_variance_list_1))/len(jackknife_variance_list)
            b = (len(RD_list_0)*np.var(RD_list_0)+len(RD_list_1)*np.var(RD_list_1))/len(RD_list)
            da = ordinary_jackknife_variance_groups(jackknife_variance_list_0,jackknife_variance_list_1)
            db = ordinary_variance_jackknife_variance_groups(RD_list_0,RD_list_1)

            uncert1,uncert2 = uncertainty(a,b,da,db)

            text_file.write(repr(a))
            text_file.write(',')
            text_file.write(repr(b))
            text_file.write(',')
            text_file.write(repr(da))
            text_file.write(',')
            text_file.write(repr(db))
            text_file.write(',')
            text_file.write(repr(uncert1))
            text_file.write(',')
            text_file.write(repr(uncert2))
            text_file.write(',')


 
        

            #Fig 4BC. Do random edits to TATACAG account for some of the variance?


            jackknife_variance_paired_list = jackknife_variance_list[np.arange(0,num_seqs,2)]+jackknife_variance_list[np.arange(0,num_seqs,2)+1]
            RD_paired_list = RD_list[np.arange(0,num_seqs,2)]-RD_list[np.arange(0,num_seqs,2)+1]

            a = np.mean(jackknife_variance_paired_list)
            b = np.var(RD_paired_list)
            da = ordinary_jackknife_variance(jackknife_variance_paired_list)
            db = ordinary_variance_jackknife_variance(RD_paired_list)

            uncert1, uncert2 = uncertainty(a,b,da,db)

            text_file.write(repr(a))
            text_file.write(',')
            text_file.write(repr(b))
            text_file.write(',')
            text_file.write(repr(da))
            text_file.write(',')
            text_file.write(repr(db))
            text_file.write(',')
            text_file.write(repr(uncert1))
            text_file.write(',')
            text_file.write(repr(uncert2))
            text_file.write(',')
    






            text_file.write('\n')
    
    
    text_file.close()

    total_time = time.clock()-pretime

if __name__ == '__main__':
     main()

