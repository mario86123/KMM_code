//
//  MixtureRankingEDA.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 04/11/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#include "MixtureRankingEDA.h"
#include <stdlib.h>
#include <limits>
#include <valarray>
#include <stdio.h>
#include <cmath>
#include "Tools.h"
#include "GeneralizedCayley.h"
#include <iomanip>

/*
 * The constructor.
 */
MixtureRankingEDA::MixtureRankingEDA(PBP * problem, int problem_size, long int max_evaluations, char * model_type, char * metric_type, int inverse, int num_clusters, char * log_filename, char * thetas_log_filename, char * weights_log_filename)
{
    //1. standard initializations
    m_problem=problem;
    m_problem_size=problem_size;
    m_max_evaluations=max_evaluations;
    m_evaluations=0;
    m_convergence_evaluations=0;

    m_pop_size=m_problem_size*10;
    m_sel_size=m_problem_size*5;
    m_offspring_size= m_pop_size-1;
    memcpy(m_metric_type,metric_type,sizeof(char)*10);
    memcpy(m_model_type,model_type,sizeof(char)*10);
    memcpy(m_log_filename,log_filename,sizeof(char)*50);
    memcpy(m_thetas_log_filename,thetas_log_filename,sizeof(char)*50);
    memcpy(m_weights_log_filename,weights_log_filename,sizeof(char)*50);
    
    m_inverse=inverse;

    //2. initialize the population
    m_population= new CPopulation(m_pop_size, m_offspring_size,m_problem_size);
    int * genes= new int[m_problem_size];
    for(int i=0;i<m_pop_size;i++)
    {
        //Create random individual
        GenerateRandomPermutation(genes,m_problem_size);
        if (m_inverse)
            m_population->SetToPopulation(genes, i, m_problem->EvaluateInv(genes));
        else
            m_population->SetToPopulation(genes, i, m_problem->Evaluate(genes));
        m_evaluations++;
    }
    delete [] genes;
    m_population->SortPopulation(0);
    m_best= new CIndividual(problem_size);
  //  m_best->SetGenes(m_population->m_individuals[0]->Genes());
  //  m_best->SetValue(m_population->m_individuals[0]->Value());

    //3. Mixtures and clusters initializations.
    m_lower_theta_bound=0.001;
    if (((string)m_metric_type)=="C")
        m_upper_theta_bound=GetUpperThetaBound_Cayley(m_problem_size);
    else if (((string)m_metric_type)=="K")
        m_upper_theta_bound=GetUpperThetaBound_Kendall(m_problem_size);

    m_EM_threshold=0.01;
    m_EM_iterations=100; //////<- usually this is 100.
    m_EM_runs=5;
    m_num_clusters=num_clusters;

    m_weights= new double[m_num_clusters];
    m_samples_models= new int[m_num_clusters];
    m_acum_weights= new double[m_num_clusters];
    m_best_weights= new double[m_num_clusters];
    
    m_z= new double*[m_num_clusters];
    for (int i=0;i<m_num_clusters;i++){
        m_z[i]= new double[m_sel_size];
    }
    m_samples= new int*[m_sel_size];
    for (int i=0;i<m_sel_size;i++){
        m_samples[i]= new int[m_problem_size];
    }
    
    m_aux_problem_size= new int[m_problem_size];
    m_aux_sel_size= new int[m_sel_size];
    for (int i=0;i<m_num_clusters;i++){
        if (((string)model_type)=="M")
            m_models.push_back(new CMallowsModel(m_problem_size, m_sel_size, metric_type, &m_lower_theta_bound, &m_upper_theta_bound));
        else if (((string)model_type)=="GM")
        {
            m_models.push_back(new CGeneralizedMallowsModel(m_problem_size, m_sel_size, metric_type, &m_lower_theta_bound, &m_upper_theta_bound));
        }
        else if (((string)model_type)=="PL")
        {
            cout<<"The implementation of Plackett-Luce mixtures is not included. "<<endl; exit(1);
        }
    }
}

/*
 * The destructor. It frees the memory allocated..
 */
MixtureRankingEDA::~MixtureRankingEDA()
{
    delete m_best;
    delete m_population;
    m_models.clear();
    
    delete [] m_weights;
    for (int i=0;i<m_num_clusters;i++)
        delete [] m_z[i];
    delete [] m_z;
    delete [] m_acum_weights;
    delete [] m_samples_models;
    delete [] m_best_weights;
    
    for (int i=0;i<m_sel_size;i++){
        delete [] m_samples[i];
    }
    delete [] m_samples;
    
    delete [] m_aux_sel_size;
    delete [] m_aux_problem_size;
}


/*
 * Returns the appropriate upper bound under the Cayley distance for the specified problem size.
 */
double MixtureRankingEDA::GetUpperThetaBound_Cayley(int problem_size){
    double upper;
    switch (problem_size) {
        case 10:
            upper=6.1;
            break;
        case 20:
            upper=7.5;
            break;
        case 30:
            upper=8.4;
            break;
        case 40:
            upper=8.9;
            break;
        case 50:
            upper=9.4;
            break;
        case 60:
            upper=9.8;
            break;
        case 70:
            upper=10.1;
            break;
        case 80:
            upper=10.3;
            break;
        case 90:
            upper=10.6;
            break;
        case 100:
            upper=10.8;
            break;
        case 150:
            upper=11.6;
            break;
        case 200:
            upper=12.2;
            break;
        case 250:
            upper=12.6;
            break;
        case 300:
            upper=13;
            break;
        case 350:
            upper=13.3;
            break;
        case 400:
            upper=13.6;
            break;
        case 450:
            upper=13.8;
            break;
        case 500:
            upper=14;
            break;
        default:
            upper=14;
            break;
    }
    return upper;
}


/*
 * Returns the appropriate upper bound under the Kendall distance for the specified problem size.
 */
double MixtureRankingEDA::GetUpperThetaBound_Kendall(int problem_size){
    double upper;
    switch (problem_size) {
        case 10:
            upper=6.1;
            break;
        case 20:
            upper=7.5;
            break;
        case 30:
            upper=8.4;
            break;
        case 40:
            upper=8.9;
            break;
        case 50:
            upper=9.4;
            break;
        case 60:
            upper=9.8;
            break;
        case 70:
            upper=10.1;
            break;
        case 80:
            upper=10.3;
            break;
        case 90:
            upper=10.6;
            break;
        case 100:
            upper=10.8;
            break;
        case 150:
            upper=11.6;
            break;
        case 200:
            upper=12.2;
            break;
        case 250:
            upper=12.6;
            break;
        case 300:
            upper=13;
            break;
        case 350:
            upper=13.3;
            break;
        case 400:
            upper=13.6;
            break;
        case 450:
            upper=13.8;
            break;
        case 500:
            upper=14;
            break;
        default:
            upper=14;
            break;
    }
    return upper;
}
/*
 * Running function
 */
int MixtureRankingEDA::Run(){
    
    //variable initializations.
    double newScore,best;
    double * weights;
    double rate=0.1;
    #ifdef PRINT_BEST
    ofstream output_file_log;
    output_file_log.open(m_log_filename);
#endif
    #ifdef PRINT_WEIGHTS
    ofstream output_file_weights;
    output_file_weights.open(m_weights_log_filename);
    double * equiprobable_weights= new double[m_num_clusters];
    for (int i=0;i<m_num_clusters;i++)
        equiprobable_weights[i]=(double)1/(double)m_num_clusters;
    #endif
    #ifdef PRINT_AVG_THETAS
    ofstream output_file_thetas;
    output_file_thetas.open(m_thetas_log_filename);
    double * averages_clusters= new double[m_num_clusters];
    #endif
    //EDA iteration. Stopping criterion is the maximum number of evaluations performed
    int iterations=0;
    while (m_evaluations<m_max_evaluations){
       
        //learn the mixture
        //cout<<"learning..."<<endl;
        weights=Learn(m_population,m_sel_size);

        //sample the mixture
        //cout<<"sampling..."<<endl;
        if ((m_evaluations+m_offspring_size)<m_max_evaluations){
            if (SAMPLING_TYPE==1)
                m_evaluations+= Sample(m_population,m_offspring_size,m_problem);
            else
                m_evaluations+= Sample_EQ(m_population,m_offspring_size,m_problem);
        }
        else{
            if (SAMPLING_TYPE==1)
                m_evaluations+= Sample(m_population,(int)(m_max_evaluations-m_evaluations),m_problem);
            else
                m_evaluations+= Sample_EQ(m_population,(int)(m_max_evaluations-m_evaluations),m_problem);
        }
                
       // cout<<"sorting..."<<endl;
        //sort the population.
        m_population->SortPopulation(1);
        

        //update indicators
        newScore=m_population->m_individuals[0]->Value();
        best=m_best->Value();
        #ifdef PRINT_AVG_THETAS
        for (int i=0;i<m_num_clusters;i++){
            averages_clusters[i]=Average(((CGeneralizedMallowsModel*)m_models[i])->m_distance_model->m_theta_parameters,m_problem_size-1);
        }
        double average_of_averages=Average(averages_clusters, m_num_clusters);
        output_file_thetas<<average_of_averages/GetUpperThetaBound(m_problem_size)<<","<<average_of_averages<<","<<Variance(averages_clusters,m_num_clusters,average_of_averages)<<endl;
        #endif
        
        #ifdef PRINT_WEIGHTS
        output_file_weights<<KullbackLeibelerDivergence(weights, equiprobable_weights, m_num_clusters)<<","<<TotalVariationDivergence(weights, equiprobable_weights, m_num_clusters)<<",";
        for (int i=0;i<m_num_clusters;i++){
            output_file_weights<<weights[i]<<",";
        }
        output_file_weights<<endl;
        #endif
        
        #ifdef PRINT_BEST
        output_file_log<<setprecision(10)<<m_population->m_individuals[0]->Value()<<endl;
        #endif
        if (newScore>best)
        {
            m_best->SetGenes(m_population->m_individuals[0]->Genes());
            m_best->SetValue(newScore);
            m_convergence_evaluations=m_evaluations;
            #ifdef PRINT_LOG
            cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<",  lower_theta: "<<m_lower_theta_bound<<", upper_theta: "<<m_upper_theta_bound<<" thetas: ";
            if (((string)m_model_type)=="M"){
                for (int i=0;i<m_num_clusters;i++)
                    cout<<((CMallowsModel*)m_models[i])->m_distance_model->m_theta_parameter<<",  ";
            }
            else{
                for (int i=0;i<m_num_clusters;i++)
                    cout<<((CGeneralizedMallowsModel*)m_models[i])->m_distance_model->m_theta_parameters[0]<<",  ";
            }
            cout<<" weights: ";
            PrintArray(weights, m_num_clusters, "");
        
            #endif
            m_lower_theta_bound=0.001;
        }
        else{
            m_lower_theta_bound+=rate;
     
            
        }
        if (m_lower_theta_bound>m_upper_theta_bound){
            m_lower_theta_bound=m_upper_theta_bound;
        }
        iterations++;
    }

   #ifdef PRINT_AVG_THETAS
    output_file_thetas.close();
    delete [] averages_clusters;
    #endif
    #ifdef PRINT_WEIGHTS
    output_file_weights.close();
    delete [] equiprobable_weights;
    #endif
    #ifdef PRINT_BEST
    output_file_log.close();
    #endif
    return 0;
}

int MixtureRankingEDA::CalculateDistanceAndX(int * sigma, int *x, int m_problem_size){
    //also updates the x vector if it isnot null
    bool * m_visited= new bool[m_problem_size];
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
    delete [] m_visited;
    return (m_problem_size-num_cycles);
}
double MixtureRankingEDA::CalculateAverageDistace_InPopulation(CPopulation * population, int sel_size){
    int i,j;
    double dist=0;
    int * m_inverted= new int[m_problem_size];
    int * m_composed= new int[m_problem_size];
    int * m_x= new int[m_problem_size-1];
    for (i=0;i<sel_size;i++){
        Invert(population->m_individuals[i]->Genes(),m_problem_size,m_inverted);
        for (j=0;j<sel_size;j++){
            Compose(population->m_individuals[j]->Genes(), m_inverted, m_composed, m_problem_size);
            dist+=CalculateDistanceAndX(m_composed, NULL, m_problem_size);
        }
    }
    delete [] m_inverted;
    delete [] m_composed;
    delete [] m_x;
    return dist/(m_sel_size);
}

/*
 * Learns the mixture model by means of the Estimation-Maximization algorithm. Returns the number of clusters generated.
 */
double * MixtureRankingEDA::Learn(CPopulation * population, int sel_size){
    
    for (int i=0;i<sel_size;i++){
        memcpy(m_samples[i], population->m_individuals[i]->Genes(),sizeof(int)*m_problem_size);
    }
    
    /* if the number of runs is 1 */
    EM(m_samples,sel_size);
    
    /* otherwise */

  /*   double loglikelihood;
     double best_loglikelihood=MIN_INTEGER;

    for (int i=0;i<m_EM_runs;i++){
        loglikelihood=EM(m_samples,sel_size);
        if (loglikelihood>best_loglikelihood){
            memcpy(m_best_weights,m_weights, sizeof(double)*m_num_clusters);
            best_loglikelihood=loglikelihood;
        }
    }

     memcpy(m_weights, m_best_weights, sizeof(double)*m_num_clusters);
*/
    return m_weights;
}


/*
 * KMeans clustering technique.
 */
void MixtureRankingEDA::KMeans(int ** samples, int sel_size, vector<CRankingModel *> models, int num_clusters, int * centroids)
{
    //Initialize arrays.
    int * centroids_indices_before= new int[num_clusters];
    
    int * belongs_to_cluster= new int[sel_size];
    int i,j,index, cluster;
    int max_iterations=10;
    vector<int> indices;
    for (i=0;i<sel_size;i++) indices.push_back(i);
    
    //1. Initializa randomly the centroids.
    for (i=0;i<num_clusters;i++){
        index=rand() % (sel_size-i);
        centroids[i]=indices.at(index);
        indices.erase(indices.begin()+index);
   //     cout<<"centroid: "<<centroids_indices_after[i]<<"   -->  "; PrintArray(samples[centroids_indices_after[i]], m_problem_size, "");
    }
    int iterations=0;
    while (memcmp(centroids_indices_before, centroids, sizeof(int)*num_clusters)!=0 && iterations <max_iterations){
        
        memcpy(centroids_indices_before, centroids, sizeof(int)*num_clusters);
        
        //2. Generate the clusters. Assign each permutation to the nearest centroid.
        int min_dist;
        int min_centroid;
        int dist, pos;
        for (i=0;i<sel_size;i++){
            min_dist=INT_MAX;
            min_centroid=0;
            if ((pos=Find(centroids_indices_before,num_clusters,i))==-1){
                for (j=0;j<num_clusters;j++){
                    dist=Cayley(samples[i], samples[centroids_indices_before[j]], m_problem_size);
                    if (dist<min_dist){
                        min_centroid=j;
                        min_dist=dist;
                    }
                }
                belongs_to_cluster[i]=min_centroid;
            }
            else{
                belongs_to_cluster[i]=pos;
            }
        }
    
//        PrintArray(belongs_to_cluster, sel_size, "belongs to: ");
    
        //3. Recalculate for each cluster the centroid.
        for (cluster=0;cluster<num_clusters;cluster++){
            centroids[cluster]=CalculateCentroid(cluster,belongs_to_cluster,samples,sel_size);
        }
//        PrintArray(centroids_indices_after, num_clusters, "centroids after: ");
//        PrintArray(centroids_indices_before, num_clusters, "centroids before: ");
        iterations++;
    }
    CPopulation * pop;
    for (i=0;i<m_num_clusters;i++){
        pop= new CPopulation(m_pop_size, 0,m_problem_size);
        index=0;
        for (j=0;j<sel_size;j++)
            if (belongs_to_cluster[j]==i){
                pop->SetToPopulation(samples[j], index, NULL);
                index++;
            }
        models[i]->Learn(pop, index);
        delete pop;
    }
    
    
    //release space.
    delete [] centroids_indices_before;
    delete [] belongs_to_cluster;
//    exit(1);
}

/*
 * Calculates the centroid of a given cluster from a set of solutions
 */
int MixtureRankingEDA::CalculateCentroid (int cluster, int * belongs_to_cluster, int **samples, int sample_size){
 //   cout<<"------------"<<endl;
 //   cout<<"cluster: "<<cluster<<endl;
    //Option 1. Choose the solution that minimizes the distance to the sample.
    int best_index;
    int best_distance=MAX_INTEGER;
    int dist;
    int i,j;
    for (i=0;i<sample_size;i++){
        if (belongs_to_cluster[i]==cluster){
            dist=0;
            for (j=0;j<sample_size;j++){
                if (belongs_to_cluster[j]==cluster){
                    dist+=Cayley(samples[i], samples[j], m_problem_size);
                }
            }
          //  cout<<"i: "<<i<<" dist: "<<dist<<endl;
            if (dist<best_distance)
            {
                best_distance=dist;
                best_index=i;
            }
        }
    }
 //   cout<<"best_distance: "<<best_distance<<" best_index: "<<best_index<<endl;
    return best_index;
}

/*
 * Executes the Estimation-Maximization algorithm and returns the loglikelihood of the estimated parameters
 */
double MixtureRankingEDA::EM(int ** samples, int sel_size)
{
    //Initialize structures.
    //0.- Variables
    int i,j,g,k;
    int ind;
    int iterations=0;
    double new_likelihood=0.0, old_likelihood=0.0, likelihood_difference=100;
    GeneralizedDistance_Model * model, *reference;
    double log_const;
double * m_ratios= new double[m_num_clusters];
    double aux_probability, acumul_probability;
    
    //1.- Weights to 1/G
    double val=(double)1/(double)m_num_clusters;
    for (i=0;i<m_num_clusters;i++)
        m_weights[i]=val;

    for (i=0;i<m_num_clusters;i++){
        
        //2.- Initialize the latent variables z to 0.
        for (j=0;j<m_sel_size;j++){
            m_z[i][j]=0;
        }
    }

    //3.- Assign the selected individuals to random clusters.
    for (i=0;i<m_sel_size;i++){
        ind= rand()% m_num_clusters;
        m_z[ind][i]=1;
    }
    
    //4.- Initialize the models.
    if (EM_INITIALIZATION==1){
        //Choose randomly the central permutations and set to 0.01 the spread parameters.
        for (i=0;i<m_num_clusters;i++)
        {
            GenerateRandomPermutation(m_aux_problem_size, m_problem_size);
            m_models[i]->Learn(m_aux_problem_size, 0.01);
        }
    }
    else{
        //K-means to calculate the central permutations.
        int * centroids = new int[m_num_clusters];
        KMeans(samples, m_sel_size, m_models, m_num_clusters, centroids);
        /*for (i=0;i<m_num_clusters;i++)
        {
            m_models[i]->Learn(samples[centroids[i]], 0.01);
        }*/
        delete [] centroids;
    }
    
    //5.- Compute initial likelihood.
    new_likelihood=0.0;
    for (i=0;i<m_sel_size;i++){
        for (g=0;g<m_num_clusters;g++){
            //new_likelihood+=( m_z[g][i] * (log(m_weights[g]/m_models[g]->NormalizationConstant()) - m_models[g]->Exponent(samples[i])) );
            log_const=0.0;
            for (j=0;j<m_problem_size-1;j++)
                log_const+=log (((CGeneralizedMallowsModel*)m_models[g])->m_distance_model->m_psis[j]);
            new_likelihood+=( m_z[g][i] * (log(m_weights[g])-log_const - m_models[g]->Exponent(samples[i])) );
        }
    }
    bool stop=false;
    //PrintMatrix(samples, m_sel_size, m_problem_size, "pop. ");
    //6.- Run the EM algorithm
    while (likelihood_difference>m_EM_threshold && iterations<m_EM_iterations ){
        
        //PrintArray(m_weights, m_num_clusters, "m_weights: ");
        //6.1- E-step. Compute the new z values.
        for (i=0;i<m_sel_size;i++){
            
            //OPTION 1: STANDARD PROCEDURE
#ifndef NUMERICAL_ERRORS
            acumul_probability=0;
            //calculate probabilities.
            for (g=0;g<m_num_clusters;g++){
                aux_probability=m_weights[g]*m_models[g]->Probability(samples[i]);
                m_z[g][i]=MAX(aux_probability,std::numeric_limits< double >::min());
                acumul_probability+=m_z[g][i];
            }
            //normalize probabilities.
            for (g=0;g<m_num_clusters;g++){
                m_z[g][i]=m_z[g][i]/acumul_probability;
            }
#else

            //OPTION 2: TO AVOID NUMERICAL ERRORS.
            reference=((CGeneralizedMallowsModel *)m_models[0])->m_distance_model;
            double ratio_acum=1.0;
            m_ratios[0]=1;
            for (g=1;g<m_num_clusters;g++){
                model=((CGeneralizedMallowsModel *)m_models[g])->m_distance_model;
                m_ratios[g]=(m_weights[g]/m_weights[0])*((GeneralizedCayley_Model*)model)->Divided_Probability(samples[i], (GeneralizedCayley_Model*)reference);
                ratio_acum+=m_ratios[g];
            }
            for (g=0;g<m_num_clusters;g++){
                m_z[g][i]=m_ratios[g]/ratio_acum;
                if (isnan(m_z[g][i])){
                    stop=true;
                }
                    
            }
#endif
        }
        if (stop==true) break;
        //6.2- M-step.
        
        //6.2.1. Compute w values.
        for (i=0;i<m_num_clusters;i++)
            m_weights[i]=0;
        for (g=0;g<m_num_clusters;g++){
            for (i=0;i<m_sel_size;i++){
                m_weights[g]+=m_z[g][i];
            }
            m_weights[g]=m_weights[g]/(double)m_sel_size;
        }

        //6.2.2. Compute the central permutations and the spread parameters taking into account the vector w of weights.
        for (k=0;k<m_sel_size;k++)
            m_aux_sel_size[k]=0;

        for (g=0;g<m_num_clusters;g++){
            m_models[g]->Learn(samples, m_sel_size, m_z[g], m_aux_sel_size);
           //m_models[g]->Learn(samples, m_sel_size, m_z[g], NULL);
        }
     
        //6.3. Compute the likelihood difference
        old_likelihood=new_likelihood;
        new_likelihood=0.0;
        for (i=0;i<m_sel_size;i++){
            for (g=0;g<m_num_clusters;g++){
                log_const=0.0;
                for (j=0;j<m_problem_size-1;j++)
                    log_const+=log (((CGeneralizedMallowsModel*)m_models[g])->m_distance_model->m_psis[j]);
                new_likelihood+=( m_z[g][i] * (log(m_weights[g])-log_const - m_models[g]->Exponent(samples[i])) );
            }
        }
     //   cout<<"new: "<<(double)new_likelihood<<"   old: "<<(double)old_likelihood<<endl;
        likelihood_difference=std::abs(new_likelihood-old_likelihood);
        iterations++;
       // PrintArray(m_weights, m_num_clusters, "");
     //   cout<<"iterations: "<<iterations<<",    "<<likelihood_difference<<endl;
    }
    delete [] m_ratios;
    //PrintArray(m_weights, m_num_clusters, "we: ");
  //  cout<<"----------"<<endl;exit(1);
    return new_likelihood;
}

/*
 * Samples solutions from the mixture model learnt using the Stochastic Universal Sampling (SUS). Returns the number of sampled solutions.
 */
int MixtureRankingEDA::Sample(CPopulation * population, int num_samples, PBP * problem)
{

    //0.- Variable initializations.
    int i,g;
    std::fill_n(m_samples_models,m_num_clusters,0);

    //SUS
    m_acum_weights[0]=m_weights[0];
    for (i=1;i<m_num_clusters;i++)
        m_acum_weights[i]=m_weights[i]+m_acum_weights[i-1];
    m_acum_weights[m_num_clusters-1]=1;// para redondear.

    //1.- Calculate with Stochastic Universal Sampling the number of individuals to sample from each model
    double ratio=(double)1/(double)num_samples;
    double acum=((double)rand()/(double)RAND_MAX+1) * ratio;

    i=0; g=0;
    while (i<num_samples){
        if (acum<=m_acum_weights[g]){
            m_samples_models[g]++;
            i++;
            acum+=ratio;
            if (acum>1) acum=1;
        }
        else{
            g++;
        }
    }
  
    //2.- Sample solutions
    for (g=0;g<m_num_clusters;g++){
        m_models[g]->Sample(population, m_samples_models[g], m_inverse, problem);
    }
    return num_samples;
}

/*
 * Samples equal number of solutions from each component in the mixture. Returns the number of sampled solutions.
 */
int MixtureRankingEDA::Sample_EQ(CPopulation * population, int num_samples, PBP * problem)
{
    
    //0.- Variable initializations.
    int i,g;
    int acum=0;
    for (i=0;i<m_num_clusters-1;i++)
    {
       m_samples_models[i]=num_samples/m_num_clusters;
        acum+=m_samples_models[i];
    }
    m_samples_models[m_num_clusters-1]=num_samples-acum;

    //2.- Sample solutions
    for (g=0;g<m_num_clusters;g++){
        m_models[g]->Sample(population, m_samples_models[g], m_inverse, problem);
    }

    return num_samples;
}

/*
 * Returns the number of performed evaluations.
 */
long int MixtureRankingEDA::GetPerformedEvaluations(){
    return m_convergence_evaluations;
}

/*
 * Returns the fitness of the best solution obtained.
 */
double MixtureRankingEDA::GetBestSolutionFitness(){
    return m_best->Value();
}

/*
 * Returns the best solution obtained.
 */
CIndividual * MixtureRankingEDA::GetBestSolution(){
    return m_best;
}























