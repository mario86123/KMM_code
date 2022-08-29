//
//  MixtureRankingEDA.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 04/11/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__MixtureRankingEDA__
#define __RankingEDAsCEC__MixtureRankingEDA__

#include <stdio.h>
#include "PBP.h"
#include "Individual.h"
#include "Population.h"
#include "RankingModel.h"
#include "MallowsModel.h"
#include "GeneralizedMallowsModel.h"
#include "Tools.h"

class MixtureRankingEDA
{
    
public:
    
    /*
     * Number of clusters.
     */
    int m_num_clusters;
    
    /*
     * Weights of clusters.
     */
    double * m_weights;
    
    /*
     * Defines the z latent variables for the EM.
     */
    double ** m_z;
    
    
    /*
     * Contains the number of solutions to sample from each model in the mixture.
     */
    int * m_samples_models;
    
    /*
     * Auxiliary array for sampling.
     */
    double * m_acum_weights;
    
    /*
     * Auxiliary array for learning.
     */
    double * m_best_weights;
    
    /*
     * Matrix of ints that represents the population.
     */
    int ** m_samples;
    
    /*
     * Auxiliary arrays for the EM.
     */
    int * m_aux_problem_size;
    int * m_aux_sel_size;

    /*
     * The size of the selection pool.
     */
    int m_sel_size;
    
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
     * The name of the file to store the logs of the weights of the clusters.
     */
    char m_weights_log_filename[50];
    
    /*
     * Evaluates the inverse of the samples solutions.
     */
    int m_inverse;
    
    /*
     * The population
     */
    CPopulation * m_population;
        /*
     * Likelihood convergence threshold of the EM.
     */
    double m_EM_threshold;

    /*
     * Maximum number of iterations inside the EM.
     */
    int m_EM_iterations;
    
    /*
     * Number of repetitions of the EM procedure.
     */
    int m_EM_runs;
    
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
     * The constructor.
     */
    MixtureRankingEDA(PBP * problem, int problem_size, long int max_evaluations, char * model_type, char * metric_type, int inverse, int num_clusters, char * log_filename, char * thetas_log_filename, char * weights_log_filename);
    
    /*
     * The destructor.
     */
    virtual ~MixtureRankingEDA();
    
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
     * Returns the appropriate upper bound under the Kendall distance for the specified problem size.
     */
    double GetUpperThetaBound_Kendall(int problem_size);
    
    /*
     * Learns the mixture model by means of the Estimation-Maximization algorithm. Returns the weights of the clusters generated.
     */
    double * Learn(CPopulation * population, int sel_size);
    
    /*
     * Executes the Estimation-Maximization algorithm and returns the loglikelihood of the estimated parameters
     */
    double EM(int ** samples, int sel_size);
    
    /*
     * Calculates the centroid of a given cluster from a set of solutions
     */
    int CalculateCentroid (int cluster, int * belongs_to_cluster, int **samples, int sample_size);
    
    /*
     * Kmeans clustering technique.
     */
    void KMeans(int ** samples, int sel_size, vector<CRankingModel *> models, int num_clusters, int * centroids);
    
    /*
     * Samples solutions from the mixture model learnt using the Stochastic Universal Sampling (SUS). Returns the number of sampled solutions.
     */
    int Sample(CPopulation * population, int num_samples, PBP * problem);
    
    /*
     * Samples equal number of solutions from each component in the mixture. Returns the number of sampled solutions.
     */
    int Sample_EQ(CPopulation * population, int num_samples, PBP * problem);
    int CalculateDistanceAndX(int * sigma, int *x, int m_problem_size);
    double CalculateAverageDistace_InPopulation(CPopulation * population, int sel_size);
};

#endif /* defined(__RankingEDAsCEC__MixtureRankingEDA__) */
