#inheritance algorithm for stochastic leakage, deterministic reamplification and homoplasmic ICs

import numpy as np
import bisect
import pandas as pd

def Roulette_wheel(I, fitnesses):
        #cell_probabilities=[f1/population_fitness, f2/population_fitness,...,fN/population_fitness] prob of being selected
        I_fathers=I[int(len(I)/2):]
        fitnesses_mothers=fitnesses[0:int(len(I)/2)]
        fitnesses_fathers=fitnesses[int(len(I)/2):]
        cum_mothers=np.cumsum(fitnesses_mothers) #cummulative sum of the fitness vector
        population_m_fitness=sum(fitnesses_mothers) #sum of all elements in the fitness vector
        X_mothers=cum_mothers/population_m_fitness
        X_mothers=X_mothers.tolist() #transform array to list
        cum_fathers=np.cumsum(fitnesses_fathers) #cummulative sum of the fitness vector
        population_f_fitness=sum(fitnesses_fathers) #sum of all elements in the fitness vector
        X_fathers=cum_fathers/population_f_fitness
        X_fathers=X_fathers.tolist() #transform array to list
        #select randomly a mother and a father
        selected_I=[]
        selected_fitnesses = []
        index_list=[]
        #mother
        r=rng.uniform(low = 0.0, high = 1, size = 1)
        f_nor_m=bisect.bisect_right(X_mothers,r) #index of the least upper bound of r in X
        if f_nor_m==len(I)/2:
            f_nor_m=f_nor_m-1
        index_list.append(f_nor_m) #index corresponds to name of individual
        selected_fitnesses.append(fitnesses[f_nor_m]) #add the fitness of that selected cell to the list 
        selected_I.append(I[f_nor_m]) #add the selected mother to the list
        #father
        r=rng.uniform(low = 0.0, high = 1, size = 1)
        f_nor_f=bisect.bisect_right(X_fathers,r) #index of the least upper bound of r in X
        if f_nor_f==len(I_fathers):
            f_nor_f=f_nor_f-1
        index_list.append(I_fathers[f_nor_f]) #index corresponds to name of individual
        selected_fitnesses.append(fitnesses[I_fathers[f_nor_f]]) #add the fitness of that selected cell to the list 
        selected_I.append(I[I_fathers[f_nor_f]]) #add the selected mother to the list
        return(selected_I, selected_fitnesses) #[mother, father] #selected_I=indices of selected mother and father in I
#env just changing once from favouring A to favouring B
def fitness(i, genetics, E, t):
    if E==0: #just one env A
        f=genetics[i][0]+0.5*genetics[i][1] 
    elif t <E:
        f=genetics[i][0]+0.5*genetics[i][1]
    elif t>=E: 
        f=0.5*genetics[i][0]+genetics[i][1]
    return f

def algorithm(I, genetics, E, t, lam, mut, seed, N): # leakage_then_mutation
    fitnesses_generations=[]
    a_generations=[]
    b_generations=[]
    m_generations=[]
    for index, t_point in enumerate(t):
        fitnesses=[]
        for i in I:
            fitnesses.append(fitness(i, genetics, E, t_point))
        fitnesses_generations.append(fitnesses)
        new_genetics=[]
        a_generations.append([])
        b_generations.append([])
        m_generations.append([])
        
        for j in I: #new individuals
            parents=Roulette_wheel(I, fitnesses)
            mother=parents[0][0]
            father=parents[0][1]
            a_pj=rng.binomial(genetics[mother][0], (1-lam)/2)+ rng.binomial(genetics[father][0], lam/2)
            b_pj=rng.binomial(genetics[mother][1], (1-lam)/2)+ rng.binomial(genetics[father][1], lam/2)
            m_pj=rng.binomial(genetics[mother][2], (1-lam)/2)+ rng.binomial(genetics[father][2], lam/2)
            
            a_mut=rng.binomial(a_pj, mut) #mutated load
            b_mut=rng.binomial(b_pj, mut) #mutated load
            a_j=a_pj-a_mut
            b_j=b_pj-b_mut
            m_j=m_pj+a_mut+b_mut
            #reamplification
            a_j=2*a_j
            b_j=2*b_j
            m_j=2*m_j
            sum_j=a_j+b_j+m_j
            
            if sum_j==0:
                sum_j=1
            a_j=round(a_j*N/sum_j)
            b_j=round(b_j*N/sum_j)
            m_j=round(m_j*N/sum_j)
            if a_j+b_j+m_j==N-1:
                if a_j<=b_j:
                    b_j=b_j+1
                elif a_j>b_j:
                    a_j=a_j+1
            elif a_j+b_j+m_j==N+1:
                if a_j<=b_j:
                    a_j=a_j-1
                elif a_j>b_j:
                    b_j=b_j-1
            #got an error b_j<0, set any type to zero if <0
            if a_j<0:
                a_j=0
            if b_j<0:
                b_j=0
            new_genetics.append([a_j, b_j, m_j])
            a_generations[index].append(a_j)
            b_generations[index].append(b_j)
            m_generations[index].append(m_j)
        genetics=new_genetics
    #save outputs as csv files
    df_a=pd.DataFrame(a_generations, index=t, columns=I, dtype=None, copy=None).T
    df_a.to_csv('a_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_'+str(seed)+'.csv', index=df_a.index.tolist())
    df_b=pd.DataFrame(b_generations, index=t, columns=I, dtype=None, copy=None).T
    df_b.to_csv('b_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_'+str(seed)+'.csv', index=df_b.index.tolist())
    df_m=pd.DataFrame(m_generations, index=t, columns=I, dtype=None, copy=None).T
    df_m.to_csv('m_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_'+str(seed)+'.csv', index=df_m.index.tolist())
    df_f=pd.DataFrame(fitnesses_generations, index=t, columns=I, dtype=None, copy=None).T
    df_f.to_csv('f_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_'+str(seed)+'.csv', index=df_f.index.tolist())
    return (df_a, df_b, df_m, df_f)

#running simulations
seeds=[666,500,200,0,1,2,3,4,5,6]
Es=[8]
lams=[0,5e-4,1e-3,2e-3,4e-3,8e-3,1.6e-2,3.2e-2,6.4e-2,1.3e-1,2.6e-1,5.1e-1]
muts=[0.01]
Ns=[10,20,40,80,160,320,640]

I=[i for i in range (0,100)]
t=[i for i in range (0,500)]

for seed in seeds:
    rng=np.random.RandomState(seed)
    for N in Ns:
        for E in Es:
            for lam in lams:
                for mut in muts:
                    print(seed, N, E, lam, mut)
                    l_1=[[N,0,0] for i in I[0:int(len(I)/4)]]
                    l_2=[[0,N,0] for i in I[int(len(I)/4):int(len(I)/2)]]
                    l_3=[[N,0,0] for i in I[int(len(I)/2):int(3*len(I)/4)]]
                    l_4=[[0,N,0] for i in I[int(3*len(I)/4):int(len(I))]]
                    genetics= [l_1, l_2, l_3, l_4]
                    genetics = sum(genetics, [])
                    algorithm(I, genetics, E, t, lam, mut, seed, N)
    

# Merge csv files for the different seeds
outputs=['a', 'b', 'm', 'f']

for N in Ns:
    for E in Es:
        for lam in lams:
            for mut in muts:
                for o in outputs:
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_666.csv'
                    df_seed666=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_500.csv'
                    df_seed500=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_200.csv'
                    df_seed200=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_0.csv'
                    df_seed0=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_1.csv'
                    df_seed1=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_2.csv'
                    df_seed2=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_3.csv'
                    df_seed3=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_4.csv'
                    df_seed4=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_5.csv'
                    df_seed5=pd.read_csv(my_csv, index_col=0)
                    my_csv = o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_seed_6.csv'
                    df_seed6=pd.read_csv(my_csv, index_col=0)
            
                    df=pd.concat([df_seed666, df_seed500])
                    df=pd.concat([df, df_seed200])
                    df=pd.concat([df, df_seed0])
                    df=pd.concat([df, df_seed1])
                    df=pd.concat([df, df_seed2])
                    df=pd.concat([df, df_seed3])
                    df=pd.concat([df, df_seed4])
                    df=pd.concat([df, df_seed5])
                    df=pd.concat([df, df_seed6])
                    df.to_csv(o+'_env_'+str(E)+'_lambda_'+str(lam)+'_mu_'+str(mut)+'_N_'+str(N)+'_allseeds.csv', index=df.index.tolist())
