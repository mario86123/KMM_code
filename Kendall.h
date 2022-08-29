//
//  Kendall.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef _KENDALL_
#define _KENDALL_
#include "Distance.h"
#include "NewtonRaphson.h"
#include "RankingEDA.h"

typedef struct kendall_data
{
    int n;
    double * Vjsmean;
    int j;
    
}kendall_data;


class Kendall_Model : public Distance_Model
{
public:

	/*
	 * The constructor.
	 */
	Kendall_Model(int problem_size, int sel_size, double * lower_theta_bound, double * upper_theta_bound);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~Kendall_Model();
	
	/*
	 * Learns the Mallows kendall model give a population of samples.
	 */
    bool Learn(int ** samples, int size);
    
    /*
     * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
	
    /*
     * Learns a Mallows model based under the Kendall distance for the given sample of individuals and the associated weights.
     */
    bool Learn(int ** samples, int size, double * weights, int * chosen);
    
    /*
     * It samples a set of new individuals.
     */
    int Sample(CPopulation * population, int num_samples, int inverse, PBP * problem);
	
	/*
	 * Calculates the probability of the individual given the probabilistic model.
	 */
	double Probability(int * individual);
	
    /*
     * Calculates the normalization constant.
     */
    double NormalizationConstant();
    
    /*
     * Calculates the exponent.
     */
    double Exponent(int * solution);
    
private:
	
    /*
     * Matrix of V probabilities.
     */
	double ** m_vprobs;
    
    /*
     * The auxiliary vector for sampling solutions.
     */
    int * m_sampling_permutation;
    
    /*
     * Frequency matrix employed when calculating the consensus ranking.
     */
    int ** m_frequency_matrix;
    
    /*
     * Auxiliary data structures for sampling stage.
     */
    int * aux;
    int * aux_v;
    int * aux_n;
    
    /*
     * Auxiliary vector used for calculating the consensus ranking.
     */
    double * consensusVector;
    
    /*
     * Auxiliary data structures for learning stage.
     */
    double * VjsMean;
	double * VjsNonMean;
	int * Vjs;
	int * invertedB;
	int * composition;
    
    /*
     * Data structure for the generic input of parameters into NewtonRaphson.
     */
    kendall_data * m_newton_d;
    
	/*
     * Calculates de consensus permutation from the given population cases.
     */
	void CalculateConsensusRanking(int** samples, int sample_size, int* consensus_ranking);
	
    /*
     * Calculates de consensus permutation from the given population cases.
     */
    void CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensus_ranking);
    
    /*
     * Calculates de consensus permutation from the given population cases.
     */
    void CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensus_ranking, int * chosen);
    
    /*
     * Calculates the sum of Kendall distances of the given solution to the sample of solutions.
     */
    double CalculateDistancetoSample(int * solution, int ** samples, int cases_num);
    
    /*
     * Calculates the sum of Kendall distances of the given solution to the sample of solutions.
     */
    double CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, double * weights);
    
    /*
     * Calculates the Kendall distance between two permutations.
     */
    int CalculateKendallDistance(int * perm_a, int * perm_b);

    /*
     * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
     */
    double CalculateThetaParameter(int*consensus_ranking, int** samples, int samples_num);
    
    /*
     * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
     */
    double CalculateThetaParameter_WeightedSamples(int*consensus_ranking, int** samples, int samples_num, double * weights);
    
    /*
     * Calculates the theta parameter that equals to 0 the function f, by brute force.
     */
    double ThetaParameter_BruteForce(kendall_data * data);
    
    /*
     * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
     */
    void CalculatePsiConstants(double theta, double* psis);
    
	/*
	 * Generates a permutation from a v vector.
	 */
	void GeneratePermuFromV(int * v, int * permu);

    /*
     * Determines if the given indiviual is a local optima for the problem and the current distance.
     */
    bool isLocalOptima(int * indiviual, PBP * problem);
	
};

#endif
