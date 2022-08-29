//
//  GeneralizedMallowsModel.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/20/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "GeneralizedMallowsModel.h"

#include "GeneralizedKendall.h"
#include "GeneralizedCayley.h"
//#include "Cayley.h"
//#include "Ulam.h"
#include "Tools.h"
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <list>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;


/*
 * Class constructor.
 */
CGeneralizedMallowsModel::CGeneralizedMallowsModel(int problem_size, int sel_size, char * metric_type, double * lower_theta_bound, double * upper_theta_bound)
{
	m_problem_size=problem_size;
    m_population_size=sel_size;

        m_population= new int*[m_population_size];

    if (((string)metric_type)=="K"){

        m_distance_model= new GeneralizedKendall_Model(problem_size, lower_theta_bound, upper_theta_bound);
    }
    else if (((string)metric_type)=="C"){

         m_distance_model= new GeneralizedCayley_Model(problem_size, lower_theta_bound, upper_theta_bound);
    }
    else{
        cout<<"There is no Generalized Mallows implementation for the metric "<<(string)metric_type<<"."<<endl;
        exit(1);
    }
    
}

/*
 * Class destructor.
 */
CGeneralizedMallowsModel::~CGeneralizedMallowsModel()
{
	delete m_distance_model;
           delete [] m_population;
}

/*
 * Virtual learning function.
 */
bool CGeneralizedMallowsModel::Learn(CPopulation * population, int size)
{
    for(int i=0;i<size;i++){
		m_population[i] = population->m_individuals[i]->Genes();
    }
    return m_distance_model->Learn(m_population, size);
}

/*
 * Virtual learning function.
 */
bool CGeneralizedMallowsModel::Learn(CPopulation * population, int size, double * weights, int * chosen)
{
    for(int i=0;i<size;i++){
        m_population[i] = population->m_individuals[i]->Genes();
    }
    return m_distance_model->Learn(m_population, size,weights,chosen);
}

/*
 * Given a population of samples, it learns a Generalized Mallows model from the data given, taking into account the probability of each sample.
 */
bool CGeneralizedMallowsModel::Learn(int ** samples, int size, double * weights, int * chosen)
{
    return m_distance_model->Learn(samples, size, weights, chosen);
}

/*
 * Builds the Generalized Mallows model with the given CR and theta parameters.
 */
bool CGeneralizedMallowsModel::Learn(int * consensus_ranking, double theta){
    return m_distance_model->Learn(consensus_ranking,theta);
}

/*
 * From the learnt model, it samples a number of individuals.
 */
int CGeneralizedMallowsModel::Sample(CPopulation * population, int num_samples, int inverse, PBP * problem)
{
    return m_distance_model->Sample(population, num_samples,inverse, problem);
}

/*
 * Calculates the probability of the solution given the probabilistic model.
 */
double CGeneralizedMallowsModel::Probability(int * solution)
{
    return m_distance_model->Probability(solution);
}

/*
 * Calculates the normalization constant.
 */
double CGeneralizedMallowsModel::NormalizationConstant()
{
    return m_distance_model->NormalizationConstant();
}
/*
 * Calculates the exponent.
 */
double CGeneralizedMallowsModel::Exponent(int * permutation)
{
    return m_distance_model->Exponent(permutation);
}
