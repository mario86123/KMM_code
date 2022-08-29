/*
 *  Variables.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

//Integer maximum value
#define MAX_INTEGER 100000000000

//Integer minimum value
#define MIN_INTEGER -10000000

//Long integer maximum value
#define MAX_LONG_INTEGER 429496729500000

//Long integer maximum value
#define MIN_LONG_INTEGER -42949672950000

//Max operation.
#define MAX(A,B) ( (A > B) ? A : B)

//Min operation.
#define MIN(A,B) ( (A < B) ? A : B)

//Boltzmann experiment
#define BOLTZMANN_ANALYSIS 0

//Print the average theta parameter at each generation
//#define PRINT_AVG_THETAS 1

//Print the fitness value of the best solution at each generation
//#define PRINT_BEST 1

//Print the weights of the clusters at each generation.
//#define PRINT_WEIGHTS 1

//Defines the sampling type of the mixture.
// 1-> Stochastic Universal Sampling (SUS)
// 2-> Equal number of samples from each model.
#define SAMPLING_TYPE 1

//To print log.
//#define PRINT_LOG 1

//Defines the initialization of the EM.
// 1-> Random initialization of the clusters.
// 2-> K-Means initialization of the clusters.
#define EM_INITIALIZATION 1

#define NUMERICAL_ERRORS 1

#define CAYLEY_SAMPLING 2