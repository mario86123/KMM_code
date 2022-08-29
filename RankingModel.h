//
//  RankingModel.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDA__RankingModel__
#define __RankingEDA__RankingModel__

#include <iostream>
#include "Population.h"
#include "Individual.h"
#include "PBP.h"


class CRankingModel
{
public:
    
	/*
	 * The constructor.
	 */
	CRankingModel();
	
	/*
	 * The destructor.
	 */
	virtual ~CRankingModel();
	
	/*
	 * Virtual learning function.
	 */
	virtual bool Learn(CPopulation * population, int size);

    /*
     * Virtual learning function.
     */
    virtual bool Learn(int ** samples, int size, double * weights, int * chosen);
    
    /*
	 * Virtual learning function.
	 */
	virtual bool Learn(int * consensus_ranking, double theta_parameter);
    
    /*
     * Virtual sampling function.
     */
    virtual int Sample(CPopulation * population, int num_samples, int inverse, PBP * problem);
    
	/*
	 * Calculates the probability of the solution given the probabilistic model.
	 */
	virtual double Probability(int * solution);
    
    /*
     * Calculates the normalization constant.
     */
    virtual double NormalizationConstant();
    
    /*
     * Calculates the exponent.
     */
    virtual double Exponent(int * permutation);
    
private:
};

#endif /* defined(__RankingEDA__RankingModel__) */
