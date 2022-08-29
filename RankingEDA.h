/*
 *  RankingEDA.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 7/11/13.
 *  Copyright 2013 University of the Basque Country. All rights reserved.
 *
 */

#ifndef _RANKINGEDA_H__
#define _RANKINGEDA_H__

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include "PBP.h"
#include "Tools.h"
#include "Population.h"
#include "RankingModel.h"

using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

class RankingEDA
{
	
public:
	
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
     * The size of the selection pool.
     */
    int m_sel_size;
    
    /*
     * The size of the offspring population.
     */
    int m_offspring_size;

    //The type of the metric to use in the model of the EDA.
    char m_metric_type[10];

    //The type of the model used in the EDA.
    char m_model_type[10];
    
    /*
     * Evaluates the inverse of the samples solutions.
     */
    int m_inverse;
    
    /*
     * The population
     */
    CPopulation * m_population;
    
    /*
     * The probalistic model.
     */
    CRankingModel * m_model;
    
    /*
     * Lower theta bound for the estimation of the theta parameters within the model.
     */
    double m_lower_theta_bound;
    
    /*
     * Upper theta bound for the estimation of the theta parameters within the model.
     */
    double m_upper_theta_bound;
    
    /*
     * The name of the file to store the logs.
     */
    char m_log_filename[50];
    
    /*
     * The name of the file to store the logs of the thetas.
     */
    char m_thetas_log_filename[50];
    
    /*
     * The constructor.
     */
	RankingEDA(PBP * problem, int problem_size, long int max_evaluations, char * model_type, char * metric_type, int inverse, char * log_filename, char * thetas_log_filename);
	
    /*
     * The destructor.
     */
    virtual ~RankingEDA();
	
    /*
     * Running function
     */
	int Run();
    
    
    int CalculateDistanceAndX(int * sigma, int *x, int m_problem_size);
    double CalculateAverageDistace_InPopulation_(CPopulation * population, int sel_size);
    
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
     * This method applies a swap of the given i,j positions in the array.
     */
    void Swap(int * array, int i, int j);
	
};
#endif
