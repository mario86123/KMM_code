//
//  MallowsModel.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "MallowsModel.h"
#include "Kendall.h"
#include "Cayley.h"
#include "Ulam2.h"
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
CMallowsModel::CMallowsModel(int problem_size, int sel_size, char * metric_type, double* lower_theta_bound, double * upper_theta_bound)
{
	m_problem_size=problem_size;
    m_population_size=sel_size;

        m_population= new int*[m_population_size];

    if (((string)metric_type)=="K")
        m_distance_model= new Kendall_Model(problem_size,sel_size, lower_theta_bound, upper_theta_bound);
    else if (((string)metric_type)=="C"){
        m_distance_model= new Cayley_Model(problem_size,sel_size, lower_theta_bound, upper_theta_bound);
    }
    else{
       //m_distance_model = new Ulam_Model(problem_size,sel_size); //implementacion con el bratelli.
        m_distance_model = new Ulam_Model2(problem_size,sel_size,lower_theta_bound, upper_theta_bound); //implementacion con ficheros de datos.
    }
}

/*
 * Class destructor.
 */
CMallowsModel::~CMallowsModel()
{
	delete m_distance_model;

    delete [] m_population;
}

/*
 * Virtual learning function.
 */
bool CMallowsModel::Learn(CPopulation * population, int size)
{
    for(int i=0;i<size;i++){
		m_population[i] = population->m_individuals[i]->Genes();
    }
    return m_distance_model->Learn(m_population, size);
}

/*
 * Given a population of samples, it learns a Mallows model from the data given, taking into account the probability of each sample.
 */
bool CMallowsModel::Learn(int ** samples, int size, double * weights, int * chosen)
{
    return m_distance_model->Learn(samples, size, weights,chosen);
}

/*
 * Builds the Mallows model with the given CR and theta parameters.
 */
bool CMallowsModel::Learn(int * consensus_ranking, double theta){
    return m_distance_model->Learn(consensus_ranking,theta);
}

/*
 * From the learnt model, it samples a number of individuals.
 */
int CMallowsModel::Sample(CPopulation * population, int num_samples, int inverse, PBP * problem)
{
    return m_distance_model->Sample(population, num_samples, inverse, problem);
}

/*
 * Calculates the probability of the solution given the probabilistic model.
 */
double CMallowsModel::Probability(int * solution)
{
    return m_distance_model->Probability(solution);
}

/*
 * Calculates the normalization constant.
 */
double CMallowsModel::NormalizationConstant()
{
    return m_distance_model->NormalizationConstant();
}
/*
 * Calculates the exponent.
 */
double CMallowsModel::Exponent(int * permutation)
{
    return m_distance_model->Exponent(permutation);
}

