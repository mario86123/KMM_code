//
//  RankingModel.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "RankingModel.h"
#include "Population.h"


/*
 * The constructor.
 */
CRankingModel::CRankingModel()
{
	
}

/*
 * The destructor.
 */
CRankingModel::~CRankingModel()
{

}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(CPopulation * population, int size)
{
    return true;
}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(int ** samples, int size, double * weights, int * chosen)
{
    return true;
}

/*
 * Virtual learning function.
 */
bool CRankingModel::Learn(int * consensus_ranking, double theta_parameter)
{
    return true;
}

/*
 * Virtual sampling function.
 */
int CRankingModel::Sample(CPopulation * population, int num_samples, int inverse, PBP * problem)
{
    return 0;
}


/*
 * Calculates the probability of the solution given the probabilistic model.
 */
double CRankingModel::Probability(int * solution)
{
    return 0;
}

/*
 * Calculates the normalization constant.
 */
double CRankingModel::NormalizationConstant(){
    return 0;
}


/*
 * Calculates the exponent.
 */
double CRankingModel::Exponent(int * permutation){
    return 0;
}