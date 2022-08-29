//
//  GeneralizedCayley.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/22/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//

#include "GeneralizedCayley.h"
#include "Tools.h"
#include "Variables.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using std::cerr;
using std::cout;
using std::endl;


/*
 * The constructor of the model.
 */
GeneralizedCayley_Model::GeneralizedCayley_Model(int problem_size, double * lower_theta_parameter, double * upper_theta_bound)
{
	m_problem_size=problem_size;
	m_psis= new double[m_problem_size-1];
	m_consensus_ranking=new int[m_problem_size];
    m_theta_parameters=new double[m_problem_size];
    m_sampling_permutation= new int[m_problem_size];
	//reservar memoria para las probabilidades de Xs.
	m_xprobs=new double[m_problem_size-1];
    
    m_lower_theta_bound=lower_theta_parameter;
    m_upper_theta_bound=upper_theta_bound;
    
    m_composed = new int[ m_problem_size ];
    m_inverted = new int[ m_problem_size ];
    m_x= new int [m_problem_size];
    m_visited = new bool[ m_problem_size ];
    m_aux= new int [m_problem_size];
    m_x_acumul = new int[ m_problem_size ];
    m_x_acumul_w = new double[ m_problem_size ];
}

/*
 * The destructor of the model.
 */
GeneralizedCayley_Model::~GeneralizedCayley_Model()
{
	delete [] m_xprobs;
	delete [] m_psis;
	delete [] m_consensus_ranking;
    delete [] m_sampling_permutation;
    delete [] m_composed;
    delete [] m_inverted;
    delete [] m_x;
    delete [] m_aux;
    delete [] m_visited;
    delete [] m_x_acumul;
    delete [] m_x_acumul_w;
}

/*
 * Learns a Generalized Mallows model based on the Cayley distance and the individuals in the model.
 */
bool GeneralizedCayley_Model::Learn(int ** samples, int size)
{
    
	//1. calculate global consensus ranking
    CalculateConsensusRanking(samples, size, m_consensus_ranking);
    
    //2. calculate theta parameters
	CalculateThetaParameters(m_consensus_ranking, samples, size, m_theta_parameters);
  
    //3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameters, m_psis);

	//4.- Calculate Xj probabilities vector.
	for(int j=0;j<m_problem_size-1;j++)
        m_xprobs[j]=(double)1/m_psis[ j ];/*(double) exp(-m_theta_parameter)/m_psis[ j ];*/

	return true;
}


/*
 * Learns a Generalized Mallows model based on the Cayley distance and the individuals in the model.
 */
bool GeneralizedCayley_Model::Learn(int * consensus_ranking, double theta)
{
    //1.- Copy the consensus ranking
    memcpy(m_consensus_ranking, consensus_ranking, sizeof(int)*m_problem_size);
    
    //2.- Copy the theta parameter.
    for (int i=0;i<m_problem_size;i++){
        m_theta_parameters[i]= theta;
    }
    //3.- Calculate psi constant.
	CalculatePsiConstants(m_theta_parameters, m_psis);
    
	//4.- Calculate Xj probabilities vector.
	for(int j=0;j<m_problem_size-1;j++)
        m_xprobs[j]=(double)1/m_psis[ j ];/*(double) exp(-m_theta_parameter)/m_psis[ j ];*/
    
    
	return true;
}

/*
 * Learns a Generalized Mallows model based under the Cayley distance for the given sample of individuals and the associated weights.
 */
bool GeneralizedCayley_Model::Learn(int ** samples, int size, double * weights, int * chosen){
    double * wei= new double[size];
    double acumul=0;
    int i;
    for (i=0;i<size;i++){
        acumul+=weights[i];
    }
    for (i=0;i<size;i++){
        wei[i]=weights[i]/acumul;
    }
    
    //1. calculate global consensus ranking
    if (chosen==NULL)
        CalculateConsensusRanking_WeightedSamples(samples, size, wei, m_consensus_ranking);
    else
        CalculateConsensusRanking_WeightedSamples(samples, size, wei, m_consensus_ranking, chosen);
    
 //   PrintArray(m_consensus_ranking, 10, "CR: ");

    //2. calculate theta parameters
    CalculateThetaParameters_WeightedSamples(m_consensus_ranking, samples, size, m_theta_parameters, wei);
    
    //3.- Calculate psi constant.
    CalculatePsiConstants(m_theta_parameters, m_psis);

    //4.- Calculate Xj probabilities vector.
    for(int j=0;j<m_problem_size-1;j++)
        m_xprobs[j]=(double)1/m_psis[ j ];/*(double) exp(-m_theta_parameter)/m_psis[ j ];*/
    
    delete [] wei;
    
    return true;
}

/*
 * Given the consensus ranking, it samples a new individual.
 */
int GeneralizedCayley_Model::Sample(CPopulation * population, int num_samples, int inverse, PBP * problem)
{
    int i, samples;
    double fitness=0;

    for (samples=0;samples<num_samples ;samples++){
        
        for(i= 0;i < m_problem_size - 1; i ++ ){
            if(((double)rand() / RAND_MAX) < m_xprobs[i])
                m_x[ i ] = 0;
            else
                m_x[ i ] =1;
            }
        m_x[m_problem_size-1] = 0;
    
        GeneratePermuFromX(m_x, m_aux);
        //Compose(m_aux, m_consensus_ranking, m_sampling_permutation, m_problem_size);
        for(i=0;i<m_problem_size;i++)
            m_sampling_permutation[i]=m_aux[m_consensus_ranking[i]];
    
    //    if (inverse)
            fitness=problem->EvaluateInv(m_sampling_permutation);
       // else
       //     fitness=problem->Evaluate(m_sampling_permutation);

        population->AddToPopulation(m_sampling_permutation, samples, fitness);
   
    }
    return num_samples;
}


/*
 * Calculates sum of the cayley distances of the solution to the sample of solutions.
 */
double GeneralizedCayley_Model::CalculateDistancetoSample(int * solution, int ** samples, int cases_num){
    
    double dist=0;
    for (int i=0;i<cases_num;i++){
        Invert(samples[i],m_problem_size,m_inverted);
        Compose(solution, m_inverted, m_composed, m_problem_size);
        dist+=CalculateDistanceAndX(m_composed, NULL);
    }
    return dist;
}

/*
 * Calculates the sum of the weighted cayley distances of the given solution to the sample of solutions
 */
double GeneralizedCayley_Model::CalculateDistancetoSample_WeightedSamples(int * solution, int ** samples, int cases_num, double * weights){
    int i,j;
    double dist=0;
    double value;

    Invert(solution,m_problem_size,m_inverted);
    for (i=0;i<cases_num;i++){
       // Invert(samples[i],m_problem_size,m_inverted);
        Compose(samples[i], m_inverted, m_composed, m_problem_size);
        CalculateDistanceAndX(m_composed, m_x);
        value=0;
        for (j = 0; j < m_problem_size-1; j++)
            value += m_x[j];
        //cout<<"value: "<<value<<endl;
        dist=dist+(value*weights[i]);
    }
    return dist;
}

/*
 * Calculates de consensus permutation from the given population cases.
 */
void GeneralizedCayley_Model::CalculateConsensusRanking(int ** samples, int cases_num, int*consensusPermutation)
{
    int option=2;
    if (option==0){
        //Option 1. Choose the best solution
        for (int i=0;i<m_problem_size;i++)
            consensusPermutation[i]=samples[0][i];
    }
    else if (option==1){
        //Option 2. Calculate the MLE consensus. Ekhine's method.
        double fitness=0;
        Approximated_MLE_Consensus(cases_num, samples, consensusPermutation, &fitness);
    }
    else{
        //Option 3. Choose the solution that minimizes the distance to the sample.
        int best_index=0;
        int best_distance=MAX_INTEGER;
        int dist;
        for (int i=0;i<cases_num;i++){
            dist=CalculateDistancetoSample(samples[i], samples, cases_num);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
        memcpy(consensusPermutation, samples[best_index],sizeof(int)*m_problem_size);
    }
}

/*
 * Calculates de set median permutation from the given population cases and the vector of weights
 */
void GeneralizedCayley_Model::CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensusPermutation)
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
    memcpy(consensusPermutation, samples[best_index],sizeof(int)*m_problem_size);
}

/*
 * Calculates de set median permutation from the given population cases and the vector of weights
 */
void GeneralizedCayley_Model::CalculateConsensusRanking_WeightedSamples(int** samples, int sample_size, double * weights, int* consensusPermutation, int * chosen)
{
    //Option 1. Choose the solution that minimizes the distance to the sample.
    int best_index=-1;
    double best_distance=MAX_INTEGER;
    double dist;
    int i;
    for (i=0;i<sample_size;i++){
        if (chosen[i]==0){
            dist=CalculateDistancetoSample_WeightedSamples(samples[i], samples, sample_size, weights);
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
    }
    if (best_index==-1){
        best_index=rand()% sample_size;
    }

    memcpy(consensusPermutation, samples[best_index],sizeof(int)*m_problem_size);
    chosen[best_index]=1;
    bool flag=false;
    for (i=best_index+1;i<sample_size;i++){
        if (memcmp(consensusPermutation, samples[i], sizeof(int)*m_problem_size)==0){
            chosen[i]=1;
            flag=true;
        }
        else if (flag==true){
            break;
        }
    }
}
/*
 * Calculates the spread theta parameter from the m_consensus_ranking and the individuals in the population.
 */
void GeneralizedCayley_Model::CalculateThetaParameters(int*consensus, int ** samples, int cases_num, double * thetas)
{
    int i,j;

    //Option 2. Calculate the theta parameter according to the spread of the solutions around the consensus ranking.
    for(i = 0 ; i < m_problem_size ; i++){
        m_x_acumul[i]=0;
        m_inverted [ consensus[ i ]] = i;
    }
    for (i = 0 ; i < cases_num; i ++ ){
        for(j = 0 ; j < m_problem_size ; j ++)
            m_composed[ j ] = samples[ i ] [ m_inverted [ j ]];
        CalculateDistanceAndX(m_composed, m_x);
        for(j = 0 ; j < m_problem_size - 1 ; j++)
          m_x_acumul[ j ] += m_x[ j ];
    }
    //PrintArray(m_x_acumul, m_problem_size-1, "x acumul: ");
    get_theta_log_likelihood(cases_num, m_x_acumul, NULL, thetas);//theta is an array of length n

    for (i=0;i<m_problem_size-1;i++){
        thetas[i]=MAX(thetas[i],*m_lower_theta_bound);
        thetas[i]=MIN(thetas[i],*m_upper_theta_bound);
    }
 
}


/*
 * Calculates the spread theta parameters from the m_consensus_ranking and the individuals in the population.
 */
void GeneralizedCayley_Model::CalculateThetaParameters_WeightedSamples(int*consensus_ranking, int** samples, int samples_num, double * thetas, double * weights)
{
    int i,j;

    //Option 2. Calculate the theta parameter according to the spread of the solutions around the consensus ranking.
    for(i = 0 ; i < m_problem_size ; i++){
        m_x_acumul_w[ i ] = 0;
        m_inverted [ consensus_ranking[ i ]] = i;
    }

    for (i = 0 ; i < samples_num; i ++ ){
        for(j = 0 ; j < m_problem_size ; j ++)
            m_composed[ j ] = samples[ i ] [ m_inverted [ j ]];
        CalculateDistanceAndX(m_composed, m_x);

        for(j = 0 ; j < m_problem_size - 1 ; j++)
            m_x_acumul_w[ j ] += (double)m_x[ j ]*weights[i];
    }

    get_theta_log_likelihood_weighted(samples_num, m_x_acumul_w, thetas);//theta is an array of length n

    for (i=0;i<m_problem_size-1;i++){
        thetas[i]=MAX(thetas[i],*m_lower_theta_bound);
        thetas[i]=MIN(thetas[i],*m_upper_theta_bound);
    }
  
}

void GeneralizedCayley_Model::get_theta_log_likelihood_weighted(int m, double *x_acumul_avg, double *theta_estim){
    
    theta_estim[m_problem_size-1] = 0;
    double xav;
    int i,j;
    for(i = 0;i < m_problem_size - 1; i ++ ){
        xav = x_acumul_avg[ i ];
        j=i+1;
      
        if(xav != 0){
      
            theta_estim[ i ] =log(m_problem_size-j)-log(xav / (1-xav));
            if(theta_estim[ i ] < 0)theta_estim[ i ] = 0;
        }
        else{
            theta_estim[ i ] = *m_upper_theta_bound;
        }
      
    }
}

/*
 * Calculates the total Psi normalization constant from the m_theta_parameter and psi-s vector.
 */
void GeneralizedCayley_Model::CalculatePsiConstants(double * thetas, double* psi)
{
    int j;
    for (int i = 0 ; i <  m_problem_size - 1; i ++ ){
        j=i+1;
        psi[ i ] = 1 + ( m_problem_size - j) *(double) exp(-thetas[ i ]);
    }
}

/*
 * Calculates the probability of the given individual in the learnt Mallows probabilistic model with the Cayley distance.
 */
double GeneralizedCayley_Model::Probability(int * individual){
    
    
    Invert(m_consensus_ranking, m_problem_size, m_inverted);
    Compose(individual, m_inverted, m_composed, m_problem_size);
    CalculateDistanceAndX(m_composed, m_x);
    
    double probability=1;
  //  double probability_aux=0;
  //  double psi=1;

    for(int i = 0 ; i < m_problem_size - 1 ; i ++ ){
        probability *= exp(-m_theta_parameters[i]*m_x[ i ])/m_psis[ i ];
//        cout<<"exp. "<<exp(-m_theta_parameters[i]*m_x[ i ])<<"   cons. "<<m_psis[ i ]<<endl;
//        psi*=m_psis[i];
//        probability_aux += -m_theta_parameters[i]*m_x[ i ] - log(m_psis[ i ]);
    }
 /*      if (probability==0){
           cout<<"prob_aux: "<<probability_aux<<endl;
            cout<<"psi: "<<psi<<endl;
            cout<<"prob_aux: "<<exp(probability_aux)<<endl;
           cout<<"¡¡Error numérico en Generalized Cayley-Probability!!"<<endl;exit(1);
       }*/
    return probability;
}

/*
 * Calculates the divided probability of a solution between the current model and another reference model.
 * Generated due to numerical errors.
 */
double GeneralizedCayley_Model::Divided_Probability(int * individual, GeneralizedCayley_Model * reference){
    Invert(m_consensus_ranking, m_problem_size, m_inverted);
    Compose(individual, m_inverted, m_composed, m_problem_size);
    CalculateDistanceAndX(m_composed, m_x);
    int j;
    double part1=1;
    double up, down;
    for (j=1;j<=m_problem_size-1; j++){
        up=(m_problem_size-j)*exp(-reference->m_theta_parameters[j-1])+1;
        down=(m_problem_size-j)*exp(-m_theta_parameters[j-1])+1;
        part1*= up/down;
    }
    
    Invert(reference->m_consensus_ranking, m_problem_size, reference->m_inverted);
    Compose(individual, reference->m_inverted, reference->m_composed, m_problem_size);
    CalculateDistanceAndX(reference->m_composed, reference->m_x);
    
    double part2=0;
    for(j = 1 ; j <= m_problem_size - 1 ; j ++ ){
        part2+= -m_theta_parameters[j-1]*m_x[j-1] + reference->m_theta_parameters[j-1]*reference->m_x[j-1];
    }
    
    double probability=part1*exp(part2);
        
    return probability;

}

/*
 * Generates a permutation from a x vector.
 */
void GeneralizedCayley_Model::GeneratePermuFromX(int *x, int *permu){
    
    if (CAYLEY_SAMPLING==1)
    {
        
    
    int distance= 0;
    int i,j;
    
    for(i = 0; i < m_problem_size ; i++)
        distance+=x[ i ];
    // cyles[][] guarda un ciclo en cada fila. El prmer numero de cada fila
    // cycles[ i ][0] es un item que tiene x_j=0 y es el mayor del ciclo
    
    // para cada elem con x_j=1 se asigma aleatoriameente a cada alguno de los ciclos
    // con la restriccion de que pertenezca a un ciclo cuyo maximo item
    //
    //
    
    //for (int i = 0 ; i < n; i ++ ) cout<< x[ i ]<<" ";cout <<" x"<<endl;
    //v2
    //int len1=distance;//num posiciones con x[ i ] ==1
    int len0=m_problem_size-distance;//num posiciones con x[ i ] ==0
    
    int *cont = new int[len0];
    int **cycles = new int*[len0];
    for (i = 0 ; i < len0; i ++ )
    {
        cont[ i ] = 0;
        cycles[ i ]  = new int[distance+1];
        for (j=0;j<distance+1;j++){
            cycles[i][j]=-1;
        }
    }
    
    int iniCont= 0;
    //en ccycles van los indices de 0..n-1
    //al construir la permu se ponen lo indices de 1..n
    for (i = 0 ; i < m_problem_size; i ++ )
        if(x[ i ] ==0){
            cycles[iniCont][cont[iniCont]++] =i;//se meten los items desde 0..n-1
            iniCont++;
        }
    
    double prob=1, partialProb= 0;
    int cycleNum;
    for(int xind=m_problem_size-1;xind>= 0;xind--){
        if(x[xind] ==1){
            cycleNum = GetMostProbCycle(xind, cycles, len0, cont);
            prob *= partialProb;
            cycles[cycleNum][cont[cycleNum]++] =xind;//se meten los items desde 0..n-1
        }
    }
    // double permsInCyle=1;
    int cycleLength;
    int *rand;
    int *cycle;
    for (i = 0 ; i < len0; i ++ ){
        cycleLength = cont[ i ];
        rand = new int[cycleLength];
        GenerateRandomPermutation(rand, cycleLength);//rand:0..n-1
        cycle=cycles[ i ];
        permu[cycle[rand[cycleLength-1]]] = cycle[rand[0]] ; //+1 <= permu[ i ] \in {1..n}
        for(j= 0;j<cycleLength-1;j++)
            permu[cycle[rand[ j ]]] = cycle[rand[j+1]] ;////+1 <=permu[ i ] \in {1..n}
        delete []cycle;
        delete [] rand;
    }
    delete [] cont;
    delete [] cycles;
}
else{
    
    int random, aux,i;
    for (  i = 0 ; i < m_problem_size ; i++) permu[ i ] = i ;
    for (  i = 0 ; i < m_problem_size - 1 ; i++)
            if ( x [ i ] == 1 ){
                random = (i + 1) + ( rand() % (m_problem_size - i - 2 + 1) ); // random \in (i + 1)+[0..n-i-2] = [i+1..n-1]
                aux = permu[  random ];
                permu[ random ] = permu[ i ];
                permu[ i ] = aux;
            }
    }
}



int GeneralizedCayley_Model::CalculateDistanceAndX(int * sigma, int *x){
    //also updates the x vector if it isnot null
    
    if(x!=NULL)for (int i = 0 ; i < m_problem_size; i ++ )x[ i ] =1;
    
    int num_cycles=0, num_visited=0, item= 0;
    
    for (int i = 0 ; i < m_problem_size; i ++ )
        m_visited[ i ] =false;
    while(num_visited<m_problem_size){
        item=num_cycles;
        while(m_visited[item])item++;
        num_cycles++;
        int maxItemInCycle= 0;
        do{
            if(item>maxItemInCycle)maxItemInCycle=item;
            m_visited[item] =true;
            num_visited++;
            item=sigma[item];
        }while(!m_visited[item]);
        if(x!=NULL)x[maxItemInCycle] = 0;
    }
    return (m_problem_size-num_cycles);
}


/*
 * Obtains the most probable cycle when building the permutation from a given vector of X.
 */
int GeneralizedCayley_Model::GetMostProbCycle(int ind, int **cycles, int len, int *leng_cycles){
    //returns the cycle index
    int i,j;
    int cont = 0;
    while(ind>cycles[cont][0]) cont++;
    //each cycles from cycle[cont] to cycle[len-1] has a particuler
    //probability of being selected as the terget
    //this prob is
    double proba_cycles[len];
    double proba_cyclesNoAcum[len];
    
    for (i = 0 ; i < len; i ++ ) proba_cycles[ i ] = 0;

    double acumul= 0;
    double multi;
    for(i=cont ; i < len; i ++ ){
        for(j=cont;j<len;j++){
            if(j==i) multi= (double)factorial(leng_cycles[ j ]);
            else  multi= (double) factorial(leng_cycles[ j ]-1);
            
            if(proba_cycles[ i ] ==0) proba_cycles[ i ] =multi;
            else proba_cycles[ i ] *= multi;
        }
        acumul += proba_cycles[ i ];
    }
    
    for (i = 0 ; i < len; i ++ ) proba_cyclesNoAcum[ i ] =(double)proba_cycles[ i ]/acumul;
    for(i=cont ; i < len; i ++ )
        if(i!=0)proba_cycles[ i ] = (double)proba_cycles[ i ]/acumul + (double)proba_cycles[i-1];
        else proba_cycles[ i ] = (double)proba_cycles[ i ]/acumul ;
    double ran = (double)rand()/(RAND_MAX);
    
    for(i=cont ; i < len; i ++ )
        if(proba_cycles[ i ]>ran){
            //*partialProb=proba_cyclesNoAcum[ i ];
            return i;
        }
    return len-1;
}

/*
 * Generates a random permutaiton with K cycles.
 */
int * GeneralizedCayley_Model::GeneratePermutation_withKcycles(int n, int k){//generateUARPermuWithKCyles
    int *sigma = new int[n];
    if(n==k) for (int i = 0 ; i < n; i ++ )sigma[ i ] =i+1;
    else
        if(k==1){
            int *cycle = new int[n];
            GenerateRandomPermutation(cycle, n);
            
            for (int i = 0 ; i < n-1; i ++ ) sigma[cycle[ i ]] = cycle[i+1]+1;
            sigma[cycle[n-1]] =cycle[0]+1;

            delete []cycle;
        }else{
            long double  ran1 =((long double)rand() / ((double )RAND_MAX +1 ) );
            int ran2=-1;
            if(ran1 < (long double)(stirling_matrix_[ n - 1 ][ k - 1 ] / stirling_matrix_[ n ][ k ])){//stirling1(n-1, k-1) / stirling1(n,k))){////n is in a cycle by itself
                /* long double ran1 = rand() % (long double)stirling1(n,k),  ran2=-1; if(ran1 < stirling1(n-1, k-1)){ */
                int *small=GeneratePermutation_withKcycles(n-1, k-1);
                for (int i = 0 ; i < n-1; i ++ )sigma[ i ] =small[ i ];
                sigma[n-1] =n;
                delete []small;
            }else{
                int *small=GeneratePermutation_withKcycles(n-1, k);
                for (int i = 0 ; i < n-1; i ++ )sigma[ i ] =small[ i ];
                ran2 = rand() % (n-1) ;//0..n-2
                sigma[n-1] =sigma[ran2];
                sigma[ran2] =n;
                delete []small;
            }
        }
    return sigma;
}



double GeneralizedCayley_Model::get_theta_log_likelihood(int m, int *x_acumul, int *x_acumul_variation, double *theta_estim){
    float likeli= 0;
    float *psi = new  float[ m_problem_size ];
    theta_estim[m_problem_size-1] = 0;
    int x_i;
    float xav;
    int i,j;
    for(i = 0;i < m_problem_size - 1; i ++ ){
        x_i = x_acumul[ i ];
        if(x_acumul_variation != NULL)
            x_i += x_acumul_variation[ i ];
        j = i+1;
        if(x_i == 0)
            x_i=1;
        if(x_i == m)
            x_i = m-1;
        xav = (float) x_i/m;
        if(xav != 0){
            theta_estim[ i ] =log(m_problem_size-j)-log(xav / (1-xav));
            if(theta_estim[ i ] < 0)theta_estim[ i ] = 0;
            psi[ i ] =1+(m_problem_size-j)*exp(-theta_estim[ i ]);
            likeli += x_i * theta_estim[ i ] + m * log(psi[ i ]);
            if ( likeli != likeli)
                likeli= 0;
        }
        else
            theta_estim[ i ] = *m_upper_theta_bound;
    }
    delete [] psi;
    return (-1 * likeli);

}


void GeneralizedCayley_Model::Approximated_MLE_Consensus( int m, int **samples, int *sigma_0, double *best_likeli){
    int **  samples_inv = new int*[ m ];
    int ** samples_copy = new int*[ m ];
    for (int i = 0 ; i < m; i++) {
        samples_copy[ i ] = new int[ m_problem_size ];
        samples_inv [ i ] = new int[ m_problem_size ];
        for(int j = 0 ; j < m_problem_size ; j ++) {
            samples_inv[ i ][ samples[ i ][ j ]] = j;
            samples_copy[ i ][ j ] = samples[ i ][ j ];
        }
    }
    approx_consensus_gmm(m, samples_copy, samples_inv, sigma_0, best_likeli);
    
    //variable_neighborhood_search(m, samples, sigma_0 , best_likeli);
    for (int i = 0 ; i < m; i++) {
        delete [] samples_inv[ i ];
        delete [] samples_copy[ i ];
    }
    delete [] samples_copy;
    delete [] samples_inv;
}



void GeneralizedCayley_Model::approx_consensus_gmm(int m, int **samples_copy, int **samples_inv, int *sigma_0, double *best_likeli){
    // pre: sigma, sigmaInv =-1
    int *freq = new int[ m_problem_size ];
    int *x_acum = new int[ m_problem_size ];
    for (int i = 0; i < m_problem_size; i++)
        x_acum[ i ] = 0;
    int * sigma_0_inv = new int[ m_problem_size ];
    
    for (int i = 0; i < m_problem_size; i++){ sigma_0[ i ] = -1;sigma_0_inv[ i ] = -1;}
    
    for (int item = 0 ; item < m_problem_size ; item ++){
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  1"<<endl;}
        for (int i = 0; i < m_problem_size; i++) freq[ i ] = 0;
        int max_freq = 0;
        int pj = -1; // 0..n-1
        //for(int s= 0;s<m;s++){  for (int i = 0 ; i < n; i ++ )cout<<samples[ s ][ i ]<<" ";cout<<"samples  trace  2"<<endl;}
        for (int s = 0;  s < m; s++){
            freq[ samples_inv[ s ][item]] ++;
            if ((freq[ samples_inv[ s ][item]]) > max_freq){
                max_freq = freq[ samples_inv[ s ][item]];
                pj = samples_inv[ s ][item] ;
            }
        }
        //for (int i = 0 ; i < n; i ++ ){for(int j= 0;j<n;j++)cout<<freq[ i ][ j ]<<" ";cout<<endl;}
        sigma_0_inv[item] = pj ;
        sigma_0 [pj] = item ;
        for (int s = 0; s < m; s++) {
            if (samples_copy[ s ][pj] != item ) {//swap
                int x = samples_copy[ s ][pj];
                int y = samples_inv[ s ][item] ;
                samples_copy[ s ][pj] = item ;
                samples_copy[ s ][y] = x;
                samples_inv [ s ][item] = pj ;
                samples_inv [ s ][x ] = y ;
                x_acum[ item ]++;
            }
        }
    }
    
    //double *theta = new double[ m_problem_size ];
    //*best_likeli = get_theta_log_likelihood(m, x_acum, NULL, theta);
   // delete [] theta;
    delete [] x_acum;
    delete [] sigma_0_inv;
    delete [] freq;
}



/*
 * Calculates the normalization constant.
 */
double GeneralizedCayley_Model::NormalizationConstant()
{   double value=1;
    //PrintArray(m_psis, m_problem_size-1, "psis: ");
    for (int i=0;i<m_problem_size-1;i++)
        value=value*m_psis[i];
    return value;
}




/*
 * Calculates the exponent.
 */
double GeneralizedCayley_Model::Exponent(int * solution){
    
    Invert(m_consensus_ranking, m_problem_size, m_inverted);
    Compose(solution, m_inverted, m_composed, m_problem_size);
    CalculateDistanceAndX(m_composed, m_x);

    double exponent=0;
    for(int i = 0 ; i < m_problem_size - 1 ; i ++ )
        exponent+=m_theta_parameters[i]*m_x[i];
        
    return exponent;
}

