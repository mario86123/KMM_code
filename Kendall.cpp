//
//  Kendall.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "Kendall.h"
#include "Variables.h"
#include "Tools.h"
#include "NewtonRaphson.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <climits>
using std::cerr;
using std::cout;
using std::endl;

/*
 * The constructor of the model.
 */
Kendall_Model::Kendall_Model(int problem_size, int sel_size, double * lower_theta_bound, double * upper_theta_bound)
{
	m_problem_size=problem_size;

    m_sampling_permutation= new int[m_problem_size];
	m_psis= new double[m_problem_size-1];
	m_consensus_ranking=new int[m_problem_size];
    
	//reservar memoria para las probabilidades de Vs.
	m_vprobs=new double*[m_problem_size-1];
	for (int i=0;i<m_problem_size-1;i++)
	{		
		m_vprobs[i]=new double[m_problem_size];
	}
    
    
    m_sample_size=sel_size;

    
    m_frequency_matrix= new int*[m_problem_size];
    for (int i=0;i<m_problem_size;i++)
        m_frequency_matrix[i]=new int[m_problem_size];
    
    m_theta_parameter=0.000001;
    m_upper_theta_bound=upper_theta_bound;
    m_lower_theta_bound=lower_theta_bound;
    m_newton_d = (kendall_data *) malloc(sizeof(kendall_data));
    m_newton_d->n=m_problem_size;
    
    //sampling stage auxiliary data structures
    aux= new int[m_problem_size];
    aux_v= new int[m_problem_size]; //last position always 0.
    aux_n= new int[m_problem_size];
    
    //learning stage auxiliary data structures
    VjsMean= new double[m_problem_size-1];
	VjsNonMean= new double[m_problem_size-1];
	Vjs = new int[m_problem_size-1];
	invertedB = new int[m_problem_size];
	composition = new int[m_problem_size];
    consensusVector=new double[m_problem_size];
    
}

/*
 * The destructor of the model.
 */
Kendall_Model::~Kendall_Model()
{
    delete [] m_consensus_ranking;
    delete [] m_sampling_permutation;
    for (int i=0;i<m_problem_size-1;i++)
	{
		delete [] m_vprobs[i];
	}
	delete [] m_vprobs;
	delete [] m_psis;
    
    for (int i=0;i<m_problem_size;i++)
        delete [] m_frequency_matrix[i];
    delete [] m_frequency_matrix;

    
    //sampling stage auxiliary data structures
    delete [] aux;
    delete [] aux_v;
    delete [] aux_n;
    
    //learning stage auxiliary data structures
	delete [] VjsMean;
	delete [] VjsNonMean;
	delete [] Vjs;
	delete [] invertedB;
	delete [] composition;
    
    	delete[] consensusVector;
}


/*
 * Calculates de consensus permutation from the given population cases.
 */
bool Kendall_Model::Learn(int ** samples, int size)
{
    
    //1.- Calculate the consensus ranking/permutation.
    CalculateConsensusRanking(samples, size, m_consensus_ranking);

    //2.- Calculate theta parameters
	m_theta_parameter= CalculateThetaParameter(m_consensus_ranking, samples, size);

//    cout<<"Theta: "<<m_theta_parameter<<endl;
	//3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameter, m_psis);
	
	//4.- Calculate Vj probabilities matrix.
	double upper, lower;
	int j,r;
	for(j=0;j<m_problem_size-1;j++)
	{
		for(r=0;r<m_problem_size-j;r++)
		{ //v(j,i) proba dql v_j tome valor i
			upper=exp((-1) * r * m_theta_parameter);
			lower=m_psis[j];
			m_vprobs[j][r] =  upper/lower;
		}
	}

    return true;
}

/*
 * Builds the Mallows model for the Kendall distance with the given CR and theta parameters.
 */
bool Kendall_Model::Learn(int * consensus_ranking, double theta)
{
    //1.- Copy the consensus ranking
    memcpy(m_consensus_ranking, consensus_ranking, sizeof(int)*m_problem_size);
 
    //2.- Copy the theta parameter.
	m_theta_parameter= theta;
    
	//3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameter, m_psis);
	
	//4.- Calculate Vj probabilities matrix.
	double upper, lower;
	int j,r;
	for(j=0;j<m_problem_size-1;j++)
	{
		for(r=0;r<m_problem_size-j;r++)
		{ //v(j,i) proba dql v_j tome valor i
			upper=exp((-1) * r * m_theta_parameter);
			lower=m_psis[j];
			m_vprobs[j][r] =  upper/lower;
		}
	}
    
    return true;
}

/*
 * Learns a Mallows model based under the Kendall distance for the given sample of individuals and the associated weights.
 */
bool Kendall_Model::Learn(int ** samples, int size, double * weights, int * chosen){
    
    double * wei= new double[size];
    double acumul=0;
    int i;
    for (i=0;i<size;i++){
        acumul+=weights[i];
    }
    for (i=0;i<size;i++){
        wei[i]=weights[i]/acumul;
    }
    
    //1.- Calculate the consensus ranking/permutation.
    if (chosen==NULL)
        CalculateConsensusRanking_WeightedSamples(samples, size, wei, m_consensus_ranking);
    else
        CalculateConsensusRanking_WeightedSamples(samples, size, wei, m_consensus_ranking, chosen);
    

    //2.- Calculate theta parameter
	m_theta_parameter= CalculateThetaParameter_WeightedSamples(m_consensus_ranking, samples, size,wei);
    
	//3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameter, m_psis);
    
  	//4.- Calculate Vj probabilities matrix.
    double upper, lower;
    int j,r;
    for(j=0;j<m_problem_size-1;j++)
    {
        for(r=0;r<m_problem_size-j;r++)
        { //v(j,i) proba dql v_j tome valor i
            upper=exp((-1) * r * m_theta_parameter);
            lower=m_psis[j];
            m_vprobs[j][r] =  upper/lower;
        }
    }
    delete [] wei;
    return true;
}

/*
 * Calculates the probability of a given individuals in the current model.
 */
double Kendall_Model::Probability(int * solution)
{
    Invert(m_consensus_ranking, m_problem_size, invertedB);
	Compose(solution,invertedB,composition,m_problem_size);
	vVector_Fast(Vjs,composition,m_problem_size,aux_v);
	
	double dist=0;
	double Psi=1.0;
    //PrintArray(m_psis, m_problem_size-1, "psis: ");
	for (int i=0;i<m_problem_size-1;i++)
	{
		dist+=Vjs[i];
		Psi=Psi*m_psis[i];
	}
    
	double probability=exp(-m_theta_parameter*dist)/Psi;

	return probability;
}



/*
 * It samples a set of new individuals.
 */
int Kendall_Model::Sample(CPopulation * population, int num_samples, int inverse, PBP * problem)
{
   // int prob[20]={10,2,6,3,12,15,0,16,17,1,5,4,9,7,18,13,14,19,8,11};
   // cout<<problem->EvaluateInv(prob)<<endl;exit(1);
    
	int i, index, limit, s;
	double randVal, acumul;

    double fitness=0;
    for (s=0;s<num_samples ;s++){
        for (i=0;i < m_problem_size;i++) aux[i]=-1;
	
        //generate samples and calcualte likelihood
        //cout << "2.1. Sample a Vj."<<endl;
	
        for(i=0;i<m_problem_size-1;i++)
        {
            //muestreo de las n-1 posiciones del vector V que definen la permutacion
            randVal=(double)rand()/((double)RAND_MAX+1);

            acumul=m_vprobs[i][0];
            index=0;
            limit=m_problem_size-1-i;
            for(;(index<limit && acumul<randVal);index++)
                acumul += m_vprobs[i][index+1];
            aux_v[i]=index;
        }

        aux_v[m_problem_size-1]=0;
        GeneratePermuFromV(aux_v,aux);
        Compose(aux,m_consensus_ranking,m_sampling_permutation,m_problem_size); //original
        fitness=0;
        if (inverse)
            fitness= problem->EvaluateInv(m_sampling_permutation);
        else
            fitness= problem->Evaluate(m_sampling_permutation);
        population->AddToPopulation(m_sampling_permutation, s, fitness);
    }

    return num_samples;
}
/*
 * Calculates de consensus permutation from the given population cases.
 */
void Kendall_Model::CalculateConsensusRanking(int** samples, int sample_size, int* consensus_ranking)
{
    int option=1;
    if (option==0){
	int i, j;
    
	//Initialize matrix to zeros.
	for (i=0;i<m_problem_size;i++)
	{
		for (j=0;j<m_problem_size;j++)
			m_frequency_matrix[i][j]=0;
	}
	
	int geneValue, genePosition;
	
	//Fill frecuency matrix reviewing individuals in the population.
	//para cada permutación muestreada:
	for (i = 0; i < sample_size; i++)
	{
		for (j=0;j<m_problem_size;j++)
		{
			geneValue=j; //original
			genePosition=samples[i][j]; //original
			
        	m_frequency_matrix[genePosition][geneValue]++;
		}
	}
	
	//cout<<"Calculate consensus vector: "<<endl;
	//Calculate consensus vector.
	int job,position;
	double positionMean;
	for (job=0;job<m_problem_size;job++)
	{
		positionMean=0;
		for (position=0;position<m_problem_size;position++)
		{
			positionMean=positionMean+position*m_frequency_matrix[position][job];
		}

		consensusVector[job]=positionMean;
	}

        RandomKeys(m_consensus_ranking, consensusVector, m_problem_size);
    }
    else{
        //Option 1. Choose the solution that minimizes the distance to the sample.
        int best_index=0;
        int best_distance=MAX_INTEGER;
        int dist;
        for (int i=0;i<sample_size;i++){
            dist=CalculateDistancetoSample(samples[i], samples, sample_size);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
        memcpy(consensus_ranking, samples[best_index],sizeof(int)*m_problem_size);
    }
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
void Kendall_Model::CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensus_ranking)
{
        //Option 1. Choose the solution that minimizes the distance to the sample.
        int best_index=0;
        double best_distance=MAX_INTEGER;
        double dist;
        for (int i=0;i<sample_size;i++){
            dist=CalculateDistancetoSample_WeightedSamples(samples[i], samples, sample_size, weights);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
        memcpy(consensus_ranking, samples[best_index],sizeof(int)*m_problem_size);
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
void Kendall_Model::CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensus_ranking, int * chosen)
{
    //Option 1. Choose the solution that minimizes the distance to the sample.
    int best_index=-1;
    double best_distance=MAX_INTEGER;
    double dist;
    for (int i=0;i<sample_size;i++){
        if (chosen[i]==0){
            dist=CalculateDistancetoSample_WeightedSamples(samples[i], samples, sample_size, weights);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
    }
//    cout<<"best index: "<<best_index<<endl;    PrintArray(chosen, m_problem_size, "chosen: ");
    if (best_index==-1){
        for (int i=0;i<sample_size;i++){
            dist=CalculateDistancetoSample_WeightedSamples(samples[i], samples, sample_size, weights);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
    }
    
    memcpy(consensus_ranking, samples[best_index],sizeof(int)*m_problem_size);
    for (int i=0;i<sample_size;i++){
        if (memcmp(consensus_ranking, samples[i], sizeof(int)*m_problem_size)==0){
            chosen[i]=1;
        }
    }
}


/*
 * Calculates the sum of Kendall distances of the given solution to the sample of solutions.
 */
double Kendall_Model::CalculateDistancetoSample(int * solution, int ** samples, int cases_num){
    
    double dist=0;
    int i,j;
    for ( i=0;i<cases_num;i++){
        //dist+=Kendall(solution, samples[i], m_problem_size, aux_v, invertedB, composition, Vjs);
        Invert(samples[i],m_problem_size,invertedB);
        Compose(solution, invertedB, composition, m_problem_size);
        vVector_Fast(Vjs,composition,m_problem_size,aux_v);
        
        for (j = 0; j < m_problem_size-1; j++)
            dist += Vjs[j];
    }
    return dist;
}


/*
 * Calculates the sum of the weighted Kendall distances of the given solution to the sample of solutions
 */
double Kendall_Model::CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, double * weights){
    

    double dist=0;
    int i,j;
    double value=0;
    Invert(solution,m_problem_size,invertedB);
    for ( i=0;i<cases_num;i++){
        Compose(samples[i], invertedB, composition, m_problem_size);
        vVector_Fast(Vjs,composition,m_problem_size,aux_v);
        value=0;
        for (j = 0; j < m_problem_size-1; j++)
            value += Vjs[j];
        
        dist=dist+(value*weights[i]);
    }

    return dist;
}
/*
 * Theta parameter estimation function.
 */
double f_Kendall(double theta, void * d)
{
    kendall_data * data= (kendall_data*)d;
	double oper1,oper2;
	double oper=( 1 /(double)( exp(theta) -1 ) );
	double results=0;
	for (int j=1;j<data->n;j++)
	{
		oper1=( (double)(data->n - j + 1) /(double)(exp((data->n - j + 1)*theta) - 1 ) );
		oper2=data->Vjsmean[j-1];
		results=results + oper -oper1 -oper2;
	}
	return results;
}

/*
 * Theta parameter estimation function derivation.
 */
double fdev_kendall(double theta, void * d)
{
    kendall_data * data= (kendall_data*)d;
	double oper1= (-1)*exp(theta) / pow(exp(theta)-1, 2);
	double oper=0;
	for (int j=1;j<data->n;j++)
	{
		oper= oper+ oper1 +  ( pow(data->n-j+1,2) * exp(theta*(data->n-j+1)) ) /(double)pow(exp(theta*(data->n-j+1))-1,2);
	}
	return oper;
}

/*
 * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
 */
double Kendall_Model::CalculateThetaParameter(int*consensus_ranking, int** samples, int samples_num)
{
	int i,j;
    for (i=0;i<m_problem_size-1;i++) {VjsNonMean[i]=0;VjsMean[i]=0;}
    
	//Calculate Vjmean vector from the population
	Invert(consensus_ranking,m_problem_size,invertedB);
    int * individua;
	for (i=0;i<samples_num;i++)
	{
		individua=samples[i];
		Compose(individua,invertedB,composition,m_problem_size);// original
		vVector_Fast(Vjs, composition, m_problem_size,aux_v);
        
		for (j=0;j<m_problem_size-1;j++)
		{
			VjsNonMean[j]=VjsNonMean[j]+Vjs[j];
		}
	}
	
	double valuee;
	for (i=0;i<m_problem_size-1;i++)
	{
		valuee=(double)VjsNonMean[i]/samples_num;
		VjsMean[i]=valuee;
	}
    
    //calculate array of thetas.
    m_newton_d->Vjsmean=VjsMean;
   // return newton(0.0001, &f_Kendall, &fdev_kendall, ((void *)m_newton_d));
    m_theta_parameter=Newton(0.0001,*m_upper_theta_bound,&f_Kendall, &fdev_kendall, ((void *)m_newton_d));
    
    m_theta_parameter=MAX(m_theta_parameter,*m_lower_theta_bound);
    m_theta_parameter=MIN(m_theta_parameter,*m_upper_theta_bound);
    return m_theta_parameter;
}

/*
 * Calculates the theta parameter that equals to 0 the function f, by brute force.
 */
double Kendall_Model::ThetaParameter_BruteForce(kendall_data * data){

    double min_result=MAX_INTEGER;
    double result=MAX_INTEGER;
    double best_theta;
    double theta=m_theta_parameter-1;
    if (theta<0) theta=0.0000001;
    
    while(result>0){
        result=f_Kendall(theta, data);
       // cout<<"t: "<<theta<<" r: "<<result<<endl;
        
        if (result<min_result){
            best_theta=theta;
            min_result=result;
        }
        theta=theta+0.00000001;
    }
    return best_theta;
}

/*
 * Calculates the spread theta parameter from the ConsensusRanking and the individuals in the population.
 */
double Kendall_Model::CalculateThetaParameter_WeightedSamples(int*consensus_ranking, int** samples, int samples_num, double * weights)
{
	int i,j;
    for (i=0;i<m_problem_size-1;i++) {VjsMean[i]=0;VjsMean[i]=0;}
    double acumul_weights=0;
    for (int i=0;i<samples_num;i++)
        acumul_weights+=weights[i];
    
	//Calculate Vjmean vector from the population
	Invert(consensus_ranking,m_problem_size,invertedB);
    int * individua;
    double weight;
	for (i=0;i<samples_num;i++)
	{
		individua=samples[i];
		Compose(individua,invertedB,composition,m_problem_size);// original
		vVector_Fast(Vjs, composition, m_problem_size,aux_v);
        weight=weights[i]/acumul_weights;
		for (j=0;j<m_problem_size-1;j++)
		{
			VjsMean[j]=VjsMean[j]+((double)Vjs[j]*weight);
		}
	}
    
    //calculate array of thetas.
    m_newton_d->Vjsmean=VjsMean;
//    m_theta_parameter=newton(m_theta_parameter,100,&f_Kendall, &fdev_kendall, ((void *)m_newton_d));
    m_theta_parameter=Newton(0.0001,*m_upper_theta_bound,&f_Kendall, &fdev_kendall, ((void *)m_newton_d));
    
    
    m_theta_parameter=MAX(m_theta_parameter,*m_lower_theta_bound);
    m_theta_parameter=MIN(m_theta_parameter,*m_upper_theta_bound);
    return m_theta_parameter;
}

/*
 * Calculates the total Psi normalization constant from the ThetaParameter and psi-s vector.
 */
void Kendall_Model::CalculatePsiConstants(double theta, double* psi)
{
	//en el param psi se dejan los valores de las ctes d normalización psi_j
	//returns psiTotal: la multiplicacion de los n-1 psi_j
	double n=m_problem_size;
	//calculate psi
	int i;
    double j;
	for(i=0;i< n - 1;i++)
	{
		j=i+1;//el indice en el paper es desde =1 .. n-1
		psi[i] = (1.0 - exp((-1.0)*(n-j+1.0)*(theta)))/(1.0 - exp((-1.0)*theta));
	}
}

/*
 * Generates a permutation from a v vector.
 */
void Kendall_Model::GeneratePermuFromV(int*v,int*permu)
{
    int val, i, index;
    for(i = 0 ; i < m_problem_size ; i ++) aux_n[i]=i;
    for(i = 0 ; i < m_problem_size - 1 ; i++){
        val = v[i];
        index = 0;
        while( !(aux_n[ index ] != -1 && val == 0))
            if(aux_n[ index ++ ] != -1)
                val --;
        permu[ i ] = index  ;
        aux_n[ index ] = -1 ;
    }
    index=0;
    while(aux_n[ index ] == -1 )index++;
    permu[ m_problem_size - 1 ] = index ;
}

/*
 * Determines if the given indiviual is a local optima for the problem and the current distance.
 */
bool Kendall_Model::isLocalOptima(int * individual, PBP * problem){
    
    int fitness= problem->Evaluate(individual);
    int * vs= new int[m_problem_size-1];
    int * neighbor_forCompose= new int[m_problem_size];
    int * neighbor= new int[m_problem_size];
    for (int i=0;i<m_problem_size-1;i++){
        for (int j=0;j<m_problem_size-1;j++)
            vs[j]=(j==i);
        
        GeneratePermuFromV(vs, neighbor_forCompose);
        
        Compose(neighbor_forCompose,individual,neighbor,m_problem_size); //original
        if (fitness<problem->Evaluate(neighbor)){
            return false;
        }
    }
    delete [] vs;
    delete [] neighbor;
    delete [] neighbor_forCompose;
    return true;
}

/*
 * Calculates the normalization constant.
 */
double Kendall_Model::NormalizationConstant()
{   double value=1;
    for (int i=0;i<m_problem_size-1;i++)
        value=value*m_psis[i];
    return value;
}

/*
 * Calculates the exponent.
 */
double Kendall_Model::Exponent(int * solution){
    
    Invert(m_consensus_ranking, m_problem_size, invertedB);
    Compose(solution,invertedB,composition,m_problem_size);
    vVector_Fast(Vjs,composition,m_problem_size,aux_v);

    double dist=0;
    for (int i=0;i<m_problem_size-1;i++)
    {
        dist+=Vjs[i];
    }
    
    double exponent=(m_theta_parameter*dist);
    return exponent;
}






