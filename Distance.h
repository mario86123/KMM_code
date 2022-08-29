//
//  Distances.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __DISTANCE_MODEL_H__
#define __DISTANCE_MODEL_H__
#include "PBP.h"
#include "Population.h"
class Distance_Model
{
public:
	
	/*
	 * Problem size.
	 */
	int m_problem_size;
	
    /*
     * Sample size.
     */
    int m_sample_size;
    
	/*
	 * Spread parameters vector.
	 */
	double m_theta_parameter;
		
	/*
	 * The consensus ranking.
	 */
	int * m_consensus_ranking;
	
	/*
	 * Psi normalization constant.
	 */
	double * m_psis;
	
    /*
     * Defines upper theta bound.
     */
    double * m_upper_theta_bound;
    
    /*
     * Defines lower theta bound.
     */
    double * m_lower_theta_bound;
    
	/*
	 * The constructor.
	 */
	Distance_Model();
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~Distance_Model();

    /*
	 * Learns a Mallows model based on the Kendall distance and the individuals in the model.
	 */
	virtual bool Learn(int ** samples, int size)=0;
    
    /*
	 * Learns a Mallows model based on the Kendall distance, the individuals in the model, and the weights associatd to the samples.
	 */
	virtual bool Learn(int ** samples, int size, double * weights, int * chosen)=0;
    
    /*
	 * Learns a Mallows model based on the Kendall distance and consensus ranking and theta parameters given.
	 */
	virtual bool Learn(int * consensus_ranking, double theta)=0;
    
    /*
     * From the learnt model, it samples a number of individuals.
     */
	virtual int Sample(CPopulation * population, int num_samples, int inverse, PBP * problem)=0;
    
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	virtual double Probability(int * individual)=0;

    /*
     * Determines if the given indiviual is a local optima for the problem and the current distance.
     */
    virtual bool isLocalOptima(int * indiviual, PBP * problem)=0;
    
    /*
     * Calculates the normalization constant.
     */
    virtual double NormalizationConstant()=0;
    
    /*
     * Calculates the exponent.
     */
    virtual double Exponent(int * permutation)=0;
private:
		

};
#endif