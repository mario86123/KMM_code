//
//  KernelRankingEDA.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 16/01/15.
//  Copyright (c) 2015 Josu Ceberio Uribe. All rights reserved.
//

#include "KernelRankingEDA.h"
#include <stdlib.h>
#include <limits>
#include <valarray>
#include <stdio.h>
#include <cmath>
#include <iomanip>

/*
 * The constructor.
 */
KernelRankingEDA::KernelRankingEDA(PBP * problem, int problem_size, long int max_evaluations, char * model_type, char * metric_type, int inverse, char * log_filename, char * thetas_log_filename, int population_size, int selection_pressure)
{
    //1. standard initializations
    m_problem=problem;
    m_problem_size=problem_size;
    m_max_evaluations=max_evaluations;
    m_evaluations=0;
    m_convergence_evaluations=0;
    

    m_pop_size=population_size;
    // printf("m_pop_size: %d\n", m_pop_size);
    // m_pop_size=m_problem_size*10;
    m_sel_size=m_pop_size/selection_pressure;
    // printf("m_sel_size: %d\n", m_sel_size);
    m_offspring_size= m_pop_size-1;
    
    m_kernel_num=m_sel_size;
    memcpy(m_metric_type,metric_type,sizeof(char)*10);
    memcpy(m_model_type,model_type,sizeof(char)*10);
    memcpy(m_log_filename,log_filename,sizeof(char)*50);
    memcpy(m_thetas_log_filename,thetas_log_filename,sizeof(char)*50);
    
    m_inverse=inverse;
    
    m_ratios= new double[m_sel_size];
    m_samples_from_each_kernel=new int[m_kernel_num];
    
    //2. initialize the population
    m_population= new CPopulation(m_pop_size, m_offspring_size,m_problem_size);
    int * genes= new int[m_problem_size];
    for(int i=0;i<m_pop_size;i++)
    {
        //Create random individual
        GenerateRandomPermutation(genes,m_problem_size);
        if (m_inverse)
            m_population->SetToPopulation(genes, i, m_problem->EvaluateInv(genes));
        else
            m_population->SetToPopulation(genes, i, m_problem->Evaluate(genes));
        m_evaluations++;
    }
    delete [] genes;
    m_population->SortPopulation(0);
    m_best= new CIndividual(problem_size);
    
    //3. Mixtures and clusters initializations.
    if (((string)m_metric_type)=="C"){
        m_upper_theta_bound=GetUpperThetaBound_Cayley(m_problem_size);
        m_lower_theta_bound=GetLowerThetaBound_Cayley(m_problem_size);
    }
    else if (((string)m_metric_type)=="K"){
        m_upper_theta_bound=GetUpperThetaBound_Kendall(m_problem_size);
        m_lower_theta_bound=GetLowerThetaBound_Kendall(m_problem_size);
    }
    m_increase_rate=0.05;
    m_theta_parameter=m_lower_theta_bound;
    for (int i=0;i<m_kernel_num;i++){
        if (((string)model_type)=="M")
            m_models.push_back(new CMallowsModel(m_problem_size, m_sel_size, metric_type, &m_lower_theta_bound, &m_upper_theta_bound));
        else if (((string)model_type)=="GM")
        {
           cout<<"The implementation of Kernels of Generalized Mallows is not included. "<<endl; exit(1);
        }
        else if (((string)model_type)=="PL")
        {
            cout<<"The implementation of Kernels of Plackett-Luce is not included. "<<endl; exit(1);
        }
    }
}

/*
 * The destructor. It frees the memory allocated..
 */
KernelRankingEDA::~KernelRankingEDA()
{
    delete m_best;
    delete m_population;
    delete [] m_ratios;
    delete [] m_samples_from_each_kernel;
    m_models.clear();
}

double KernelRankingEDA::GetLowerThetaBound_Cayley(int problem_size){
    double upper;
    switch (problem_size) {
        case 10:
            upper=4.1;
            break;
        case 20:
            upper=5.5;
            break;
        case 30:
            upper=6.4;
            break;
        case 40:
            upper=7;
            break;
        case 50:
            upper=7.6;
            break;
        case 60:
            upper=8.2;
            break;
        case 70:
            upper=8.5;
            break;
        case 80:
            upper=8.6;
            break;
        case 90:
            upper=8.7;
            break;
        case 100:
            upper=8.8;
            break;
        case 150:
            upper=9.5;
            break;
        case 200:
            upper=10.2;
            break;
        case 250:
            upper=10.6;
            break;
        case 300:
            upper=11.1;
            break;
        case 350:
            upper=11.3;
            break;
        case 400:
            upper=11.5;
            break;
        case 450:
            upper=11.9;
            break;
        case 500:
            upper=12.1;
            break;
        default:
            upper=3;
            break;
    }
    return upper;
}
double KernelRankingEDA::GetUpperThetaBound_Cayley(int problem_size){
    double upper;
    switch (problem_size) {
        case 10:
            upper=6.1;
            break;
        case 20:
            upper=7.5;
            break;
        case 30:
            upper=8.4;
            break;
        case 40:
            upper=8.9;
            break;
        case 50:
            upper=9.4;
            break;
        case 60:
            upper=9.8;
            break;
        case 70:
            upper=10.1;
            break;
        case 80:
            upper=10.3;
            break;
        case 90:
            upper=10.6;
            break;
        case 100:
            upper=10.8;
            break;
        case 150:
            upper=11.6;
            break;
        case 200:
            upper=12.2;
            break;
        case 250:
            upper=12.6;
            break;
        case 300:
            upper=13;
            break;
        case 350:
            upper=13.3;
            break;
        case 400:
            upper=13.6;
            break;
        case 450:
            upper=13.8;
            break;
        case 500:
            upper=14;
            break;
        default:
            upper=14;
            break;
    }
    return upper;
}
double KernelRankingEDA::GetLowerThetaBound_Kendall(int problem_size){
    double upper;
    switch (problem_size) {
        case 10:
            upper=2.6;
            break;
        case 20:
            upper=3.3;
            break;
        case 30:
            upper=3.7;
            break;
        case 40:
            upper=4.1;
            break;
        case 50:
            upper=4.3;
            break;
        case 60:
            upper=4.5;
            break;
        case 70:
            upper=4.6;
            break;
        case 80:
            upper=4.7;
            break;
        case 90:
            upper=4.9;
            break;
        case 100:
            upper=5.0;
            break;
        case 150:
            upper=5.4;
            break;
        case 200:
            upper=5.6;
            break;
        case 250:
            upper=5.8;
            break;
        case 300:
            upper=6.2;
            break;
        case 350:
            upper=6.3;
            break;
        case 400:
            upper=6.4;
            break;
        case 450:
            upper=6.5;
            break;
        case 500:
            upper=6.6;
            break;
        default:
            upper=3;
            break;
    }
    return upper;
}

/*
 * Returns the appropriate upper bound under the KEndall distance for the specified problem size.
 */
double KernelRankingEDA::GetUpperThetaBound_Kendall(int problem_size){
    double upper;
    switch (problem_size) {
            case 10:
            upper=4.4;
            break;
            case 20:
            upper=5.2;
            break;
            case 30:
            upper=5.5;
            break;
            case 40:
            upper=5.8;
            break;
            case 50:
            upper=6.1;
            break;
            case 60:
            upper=6.2;
            break;
            case 70:
            upper=6.3;
            break;
            case 80:
            upper=6.4;
            break;
            case 90:
            upper=6.6;
            break;
            case 100:
            upper=6.7;
            break;
            case 150:
            upper=7.2;
            break;
            case 200:
            upper=7.5;
            break;
            case 250:
            upper=7.7;
            break;
            case 300:
            upper=8.0;
            break;
            case 350:
            upper=8.1;
            break;
            case 400:
            upper=8.2;
            break;
            case 450:
            upper=8.3;
            break;
            case 500:
            upper=8.4;
            break;
        default:
            upper=14;
            break;
    }
    return upper;
}
/*
 * Running function
 */
int KernelRankingEDA::Run(){
    
    //variable initializations.
    double newScore,best;
#ifdef PRINT_BEST
    ofstream output_file_log;
    output_file_log.open(m_log_filename);
#endif
  //  int non_improvement_iterations=0;
    //EDA iteration. Stopping criterion is the maximum number of evaluations performed
    int iterations=0;
    while (m_evaluations<m_max_evaluations){
        
        //learn the mixture
        //cout<<"learning..."<<endl;
        Learn(m_population,m_sel_size);
        
        //sample the mixture
        //cout<<"sampling..."<<endl;
        if ((m_evaluations+m_offspring_size)<m_max_evaluations){
            m_evaluations+= Sample(m_population,m_offspring_size,m_problem);
        }
        else{
            m_evaluations+= Sample(m_population,(int)(m_max_evaluations-m_evaluations),m_problem);
        }
        
        // cout<<"sorting..."<<endl;
        //sort the population.
        m_population->SortPopulation(1);
        
        //update indicators
        newScore=m_population->m_individuals[0]->Value();
        best=m_best->Value();
        
#ifdef PRINT_BEST
        output_file_log<<setprecision(10)<<m_population->m_individuals[0]->Value()<<endl;
#endif
        if (newScore>best)
        {
            m_best->SetGenes(m_population->m_individuals[0]->Genes());
            m_best->SetValue(newScore);
            m_convergence_evaluations=m_evaluations;
#ifdef PRINT_LOG
            cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<", theta: "<<m_theta_parameter<<",  lower_theta: "<<m_lower_theta_bound<<", upper_theta: "<<m_upper_theta_bound<<endl;
#endif
            m_theta_parameter=m_lower_theta_bound;
        //    non_improvement_iterations=0;
        }
        else{
            m_theta_parameter+=m_increase_rate;
          /*  non_improvement_iterations++;
            if (non_improvement_iterations>100){
                cout<<"restart"<<endl;
                Restart(m_population);
                
            }*/
        }
        if (m_theta_parameter>m_upper_theta_bound){
            m_theta_parameter=m_upper_theta_bound;
        }
        iterations++;
    }
    
#ifdef PRINT_BEST
    output_file_log.close();
#endif
    return 0;
}
void KernelRankingEDA::Restart(CPopulation * population){
    int * genes= new int[m_problem_size];
    for (int i=1;i<m_pop_size;i++){
        GenerateRandomPermutation(genes,m_problem_size);
    if (m_inverse)
        m_population->SetToPopulation(genes, i, m_problem->EvaluateInv(genes));
    else
        m_population->SetToPopulation(genes, i, m_problem->Evaluate(genes));
    m_evaluations++;
    }
    delete [] genes;
    
}
/*
 * Learns the kernels of the specified model by means.
 */
void KernelRankingEDA::Learn(CPopulation * population, int sel_size){
    
    for (int i=0;i<sel_size;i++){
        m_models[i]->Learn(population->m_individuals[i]->Genes(), m_theta_parameter);
    }
}


/*
 * Samples solutions from the mixture model learnt using the Stochastic Universal Sampling (SUS). Returns the number of sampled solutions.
 */
int KernelRankingEDA::Sample(CPopulation * population, int num_samples, PBP * problem)
{
    int i,g;
    double avg_fitness=0;
    
    for (i=0;i<m_sel_size;i++){
        avg_fitness+= INT_MAX-m_population->m_individuals[i]->Value();
        m_ratios[i]= INT_MAX-m_population->m_individuals[i]->Value();
    }
    
    for (i=0;i<m_sel_size;i++){
        m_ratios[i]=m_ratios[i]/avg_fitness;
    }
    

    int toRound=0;
    for (i=1;i<m_kernel_num;i++){
        m_samples_from_each_kernel[i]=(int)(m_ratios[i]*m_offspring_size);
        toRound+=m_samples_from_each_kernel[i];
        //m_samples_from_each_kernel[i]=1;
    }
    m_samples_from_each_kernel[0]=m_offspring_size-toRound;
    //PrintArray(m_samples_from_each_kernel, m_kernel_num, "samples: "); exit(1);
    //2.- Sample solutions
    for (g=0;g<m_kernel_num;g++){
        m_models[g]->Sample(population, m_samples_from_each_kernel[g], m_inverse, problem);
    }

    return num_samples;
}

/*
 * Returns the number of performed evaluations.
 */
long int KernelRankingEDA::GetPerformedEvaluations(){
    return m_convergence_evaluations;
}

/*
 * Returns the fitness of the best solution obtained.
 */
double KernelRankingEDA::GetBestSolutionFitness(){
    return m_best->Value();
}

/*
 * Returns the best solution obtained.
 */
CIndividual * KernelRankingEDA::GetBestSolution(){
    return m_best;
}









