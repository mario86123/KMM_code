//
//  ExperimentBoltzmann.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 04/11/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#include "ExperimentBoltzmann.h"
#include "MallowsModel.h"
#include "PlackettLuce.h"
#include "Ulam2.h"
#include "Tools.h"
#include "LOP.h"
/*
 * The constructor.
 */
CExperimentBoltzmann::CExperimentBoltzmann(){
    //This experiment first calculates the exact Boltzmann distribution of an instance of the LOP by brute force. Then,
    // a Mallows-Ulam and a Plackett-Luce are estimated from this Boltzmann distribution. Finally, the Kullback
    //Leibler divergence is calculated with respect to the Boltzmann distribution.
}

/*
 * The destructor.
 */
CExperimentBoltzmann::~CExperimentBoltzmann(){
    
}

/*
 * Running function.
 */
void CExperimentBoltzmann::Run(int problem_size, PBP * problem){
    
    //some initializations
    double m_lower_theta_bound=0.001;
    double m_upper_theta_bound=10;
    
    //Generate the group of all permutations of size n.
    int num_perms= factorial(problem_size);
    int ** search_space= new int*[num_perms];
    int p1[problem_size];
    for (int i=0;i<problem_size;i++) p1[i]=i;
    int index=0;
    sort(p1,p1+problem_size);
    do{
        search_space[index]=new int[problem_size];
        memcpy(search_space[index],p1,sizeof(int)*problem_size);
        index++;
    }
    while ( next_permutation (p1,p1+problem_size) );
    
    //Compare the distributions.
    CRankingModel * mallows_ulam= new CMallowsModel(problem_size,1,"U",&m_lower_theta_bound, &m_upper_theta_bound);
    Distance_Model * ulam=((CMallowsModel*)mallows_ulam)->m_distance_model;
    CRankingModel * plackett_luce= new CPlackettLuceModel(problem_size,num_perms);
    for (double c=10;c<=10;c=c+10){
        
        double * boltzmann_distribution=((LOP*)problem)->Calculate_BoltzmannDistribution(c);
        double * mallows_distribution= ((Ulam_Model2*)ulam)->Learn_fromBoltzmann(search_space, boltzmann_distribution,num_perms);
        double * plackett_distribution =((CPlackettLuceModel*)plackett_luce)->Learn_fromBoltzmann(search_space,boltzmann_distribution,num_perms);
        
        cout<<c<<","<<KullbackLeibelerDivergence(boltzmann_distribution, mallows_distribution, num_perms)<<","<<KullbackLeibelerDivergence(boltzmann_distribution, plackett_distribution, num_perms)<<endl;
        
        delete [] boltzmann_distribution;
        delete [] mallows_distribution;
        delete [] plackett_distribution;
    }
    delete mallows_ulam;
    delete ulam;
    delete plackett_luce;
    for (int i=0;i<num_perms;i++)
        delete [] search_space[i];
    delete [] search_space;
}