//
//  GeneralizedKendall.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/20/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__GeneralizedKendall__
#define __RankingEDAsCEC__GeneralizedKendall__

#include "GeneralizedDistance.h"


typedef struct kendall_data
{
    int n;
    double * Vjsmean;
    int j;
    
}kendall_data;


class GeneralizedKendall_Model : public GeneralizedDistance_Model
{
public:
    
	/*
	 * The constructor.
	 */
	GeneralizedKendall_Model(int problem_size, double * lower_theta_parameter, double * upper_theta_bound);
	
	/*
	 * The destructor. It frees the memory allocated at the construction of the kendall model.
	 */
	virtual ~GeneralizedKendall_Model();
	
	/*
	 * Learns a Generalized Mallows kendall model give a population of samples.
	 */
    bool Learn(int ** samples, int size);
	
    /*
     * Builds a Generalized Mallows model for the Kendall distance with the given CR and theta parameters.
     */
    bool Learn(int * consensus_ranking, double theta);
    
    /*
     * Learns a Generalized Mallows model based under the Kendall distance for the given sample of individuals and the associated weights.
     */
    bool Learn(int ** samples, int size, double * weights, int * chosen);
    
	/*
	 * Given the consensus ranking, it samples a new individual.
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
     * Frequency matrix employed when calculating the consensus ranking.
     */
    int ** m_frequency_matrix;
    
    /*
     * The auxiliary vector for sampling solutions.
     */
    int * m_sampling_permutation;
    
    /*
     * Auxiliary data structures for sampling stage.
     */
    int * aux;
    int * aux_v;
    int * aux_n;
    
    /*
     *
     */
    double * consensusVector;
    
    /*
     * Auxiliary data structures for learning stage.
     */
    double * VjsMean;
	int * VjsNonMean;
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
     * Calculates de set median permutation from the given population cases and the vector of weights
     */
    void CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensus_ranking);
    /*
     * Calculates de set median permutation from the given population cases and the vector of weights
     */
    void CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensus_ranking, int * chosen);
    
    /*
     * Calculates the sum of the weighted Kendall distances of the given solution to the sample of solutions
     */
    double CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, double * weights);
    
    /*
     * Calculates the distances of the given solution to the sample.
     */
    double CalculateDistancetoSample(int * solution, int ** samples, int cases_num);

    /*
     * Calculates the spread theta parameters from the ConsensusRanking and the individuals in the population.
     */
    void CalculateThetaParameters(int*consensus_ranking, int** samples, int samples_num,double * theta_parameters);
    
    /*
     * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population taking into account the vector of weights.
     */
    void CalculateThetaParameters_WeightedSamples(int*consensus_ranking, int** samples, int samples_num, double * weights);
    
    /*
     * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
     */
    void CalculatePsiConstants(double* thetas, double* psi);
    
	/*
	 * Generates a permutation from a v vector.
	 */
	void GeneratePermuFromV(int * v, int * permu);
    	
};

#endif /* defined(__RankingEDAsCEC__GeneralizedKendall__) */
