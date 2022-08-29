//
//  ExperimentBoltzmann2.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 10/11/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#include "ExperimentBoltzmann2.h"
#include "LOP.h"
/*
 * The constructor.
 */
CExperimentBoltzmann2::CExperimentBoltzmann2(){
    //This experiment performs the sampling of the Boltzmann distribution under the LOP, based on its L-decomposable property.
    
}

/*
 * The destructor.
 */
CExperimentBoltzmann2::~CExperimentBoltzmann2(){
    
}

/*
 * Running function.
 */
void CExperimentBoltzmann2::Run(int problem_size, LOP * problem){
   
    //PrintMatrix(problem->m_matrix, problem_size, problem_size, "LOP: ");
    double beta=0.1;
    int i,j;
    
    //initialize to 0 the marginals matrix.
    double ** marginals_matrix= new double*[problem_size];
    double ** contributions_matrix= new double*[problem_size];
    for (i=0;i<problem_size;i++){
        marginals_matrix[i]= new double[problem_size];
        contributions_matrix[i]= new double[problem_size];
        for (j=0;j<problem_size;j++){
            marginals_matrix[i][j]=0;
            contributions_matrix[i][j]=0;
        }
    }
    
   
    double aux_value=0;
    int * solution= new int[problem_size];
    for (i=0;i<problem_size;i++){
        solution[i]=i;
    }
    sort(solution,solution+problem_size);
    
    double eval=0;

    do{
        eval=(double)problem->Evaluate(solution)*beta;
        aux_value = exp(eval);
        for (i=0;i<problem_size;i++){
            marginals_matrix[solution[i]][i]+=aux_value;
            contributions_matrix[solution[i]][i]+=problem->Contribution(solution, solution[i], i);
        }
    }
    while ( next_permutation (solution,solution+problem_size) );
    delete [] solution;

  
    //normalize matrix.
    double acum_marginals=0;
    double acum_contributions=0;
    for (i=0;i<problem_size;i++){
        acum_marginals=0;
        acum_contributions=0;
        for (j=0;j<problem_size;j++){
            acum_marginals+=marginals_matrix[i][j];
            acum_contributions+=contributions_matrix[i][j];
        }
        for (j=0;j<problem_size;j++){
            marginals_matrix[i][j]=marginals_matrix[i][j]/acum_marginals;
            contributions_matrix[i][j]=contributions_matrix[i][j]/acum_contributions;
        }
    }
        PrintMatrix(contributions_matrix, problem_size, problem_size, "Contributions: ");
  PrintMatrix(marginals_matrix, problem_size, problem_size, "Marginals: ");
    
    exit(1);
}