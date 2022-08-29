//
//  KernelRankingEDA.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 16/01/15.
//  Copyright (c) 2015 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__KernelRankingEDA__
#define __RankingEDAsCEC__KernelRankingEDA__


#include <stdio.h>
#include "PBP.h"
#include "Individual.h"
#include "Population.h"
#include "RankingModel.h"
#include "MallowsModel.h"
#include "Tools.h"

class KernelRankingEDA
{
    
public:
    
    /*
     * The size of the selection pool.
     */
    int m_sel_size;

    /*
     * The number of kernels to learn.
     */
    int m_kernel_num;
    
    /*
     * Problem
     */
    PBP * m_problem;
    
    /*
     * Problem size
     */
    int m_problem_size;
    
    /*
     * Maximum number of evaluations permitted in the execution.
     */
    long int m_max_evaluations;
    
    /*
     * The number of evaluations performed at the moment.
     */
    long int m_evaluations;
    
    /*
     * The number of the evaluations performed to achieve the best solution so far.
     */
    long int m_convergence_evaluations;
    
    /*
     * The best solution found so far.
     */
    CIndividual * m_best;
    
    /*
     * The size of the population.
     */
    int m_pop_size;
    
    /*
     * The size of the offspring population.
     */
    int m_offspring_size;
    
    /*
     * The type of the metric to use in the model of the EDA.
     */
    char m_metric_type[10];
    
    /*
     * The type of the model used in the EDA.
     */
    char m_model_type[10];
    
    /*
     * The name of the file to store the logs of the thetas.
     */
    char m_log_filename[50];
    
    /*
     * The name of the file to store the logs of the thetas.
     */
    char m_thetas_log_filename[50];
    
    
    /*
     * Evaluates the inverse of the samples solutions.
     */
    int m_inverse;
    
    /*
     * Array to store the weight of each kernel to be used in the sampling.
     */
    double * m_ratios;
    
    /*
     * Array to store the number of individuals that are samples from each sample.
     */
    int * m_samples_from_each_kernel;
    
    /*
     * The population
     */
    CPopulation * m_population;

    /*
     * The probalistic model.
     */
    vector<CRankingModel *> m_models;
    
    /*
     * Lower theta bound for the estimation of the theta parameters within the model.
     */
    double m_lower_theta_bound;
    
    /*
     * Upper theta bound for the estimation of the theta parameters within the model.
     */
    double m_upper_theta_bound;

    /*
     * The theta parameter
     */
    double m_theta_parameter;
    
    /*
     * Increase rate of the temperature.
     */
    double m_increase_rate;
    /*
     * The constructor.
     */
    KernelRankingEDA(PBP * problem, int problem_size, long int max_evaluations, char * model_type, char * metric_type, int inverse, char * log_filename, char * thetas_log_filename);
    
    /*
     * The destructor.
     */
    virtual ~KernelRankingEDA();
    
    void Restart(CPopulation * population);
        
    /*
     * Running function
     */
    int Run();
    
    /*
     * Returns the number of performed evaluations.
     */
    long int GetPerformedEvaluations();
    
    /*
     * Returns the fitness of the best solution obtained.
     */
    double GetBestSolutionFitness();
    
    /*
     * Returns the best solution obtained.
     */
    CIndividual * GetBestSolution();
    
private:
    
    /*
     * Returns the appropriate upper bound under the Cayley distance for the specified problem size.
     */
    double GetUpperThetaBound_Cayley(int problem_size);
    /*
     * Returns the appropriate lower bound under the Cayley distance for the specified problem size.
     */
    double GetLowerThetaBound_Cayley(int problem_size);

    /*
     * Returns the appropriate upper bound under the Cayley distance for the specified problem size.
     */
    double GetUpperThetaBound_Kendall(int problem_size);
    /*
     * Returns the appropriate lower bound under the Cayley distance for the specified problem size.
     */
    double GetLowerThetaBound_Kendall(int problem_size);
    
    /*
     * Learns the kernels of the specified model by means.
     */
    void Learn(CPopulation * population, int sel_size);
    
    /*
     * Samples solutions from the mixture model learnt using the Stochastic Universal Sampling (SUS). Returns the number of sampled solutions.
     */
    int Sample(CPopulation * population, int num_samples, PBP * problem);
    
};

#endif /* defined(__RankingEDAsCEC__KernelRankingEDA__) */
