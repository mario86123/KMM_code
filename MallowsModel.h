//
//  MallowsModel.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDA__MallowsModel__
#define __RankingEDA__MallowsModel__

#include <iostream>

#include "Individual.h"
#include "RankingModel.h"
#include "Population.h"
#include "Distance.h"
#include <list>
#include <vector>


class CMallowsModel : public CRankingModel
{
	
public:
	
    /*
     * Distance model.
     */
    Distance_Model * m_distance_model;
    
    /*
     * The constructor.
     */
	CMallowsModel(int problem_size, int sel_size, char * metric_type, double * lower_theta_bound, double * upper_theta_bound);
	
    /*
     * The destructor.
     */
    virtual ~CMallowsModel();
	
    /*
     * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
    
    /*
     * Given a population of samples, it learns a Mallows model from the data given, taking into account the probability of each sample.
     */
    bool Learn(int ** samples, int size, double * weights, int * chosen);

    /*
     * Given a population of samples, it learns a Mallows model from the data.
     */
    bool Learn(CPopulation * population, int size);
    
    /*
     * Creates a new individual sampling the generalized mallows probabilistic model.
     */
    CIndividual * Simulate();
    
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	double Probability(int * individual);
    
    /*
     * From the learnt model, it samples a number of individuals.
     */
    int Sample(CPopulation * population, int num_samples, int inverse, PBP * problem);
	
    /*
     * Calculates the normalization constant.
     */
    double NormalizationConstant();
    
    /*
     * Calculates the exponent.
     */
    double Exponent(int * permutation);
    
private:
    
    /*
     * Problem size.
     */
    int m_problem_size;

    /*
     * Sample size.
     */
    int m_population_size;
    
    /*
     * Sample of solutions.
     */
    int ** m_population;
    
    
	
	};

#endif /* defined(__RankingEDA__MallowsModel__) */
