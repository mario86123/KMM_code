/*
 *  RankingEDA.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 9/15/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "RankingEDA.h"
#include "MallowsModel.h"
#include "Cayley.h"
#include "Kendall.h"
#include "GeneralizedMallowsModel.h"
#include "PlackettLuce.h"
#include <iomanip>
/*
 * The constructor.
 */
RankingEDA::RankingEDA(PBP * problem, int problem_size, long int max_evaluations, char* model_type, char* metric_type, int inverse,char * log_filename, char * thetas_log_filename)
{
    //1. standard initializations
    m_problem=problem;
    m_problem_size=problem_size;
    m_max_evaluations=max_evaluations;
    m_evaluations=0;
    m_convergence_evaluations=0;
    m_best= new CIndividual(problem_size);
    m_best->SetValue(MIN_LONG_INTEGER);
    m_pop_size=m_problem_size*10;
    m_sel_size=m_problem_size;
    m_offspring_size= m_pop_size-1;
    memcpy(m_metric_type,metric_type,sizeof(char)*10);
    memcpy(m_model_type,model_type,sizeof(char)*10);
    memcpy(m_thetas_log_filename,thetas_log_filename,sizeof(char)*50);
    memcpy(m_log_filename,log_filename,sizeof(char)*50);
    m_inverse=inverse;
    m_lower_theta_bound=0.001;
    if (((string)m_metric_type)=="C")
        m_upper_theta_bound=GetUpperThetaBound_Cayley(m_problem_size);
    else if (((string)m_metric_type)=="K")
        m_upper_theta_bound=GetUpperThetaBound_Kendall(m_problem_size);
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
    
    m_population->SortPopulation(0);
   // cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<endl;
    delete [] genes;
    
    
    //3. Build model structures
    if (((string)model_type)=="M"){
        m_model=new CMallowsModel(m_problem_size, m_sel_size, metric_type, &m_lower_theta_bound, &m_upper_theta_bound);
    }
    else if (((string)model_type)=="GM")
    {
        m_model=new CGeneralizedMallowsModel(m_problem_size, m_sel_size, metric_type, &m_lower_theta_bound, &m_upper_theta_bound);
    }
    else if (((string)model_type)=="PL")
    {
        m_model= new CPlackettLuceModel(m_problem_size, m_sel_size);
    }
}

/*
 * The destructor. It frees the memory allocated..
 */
RankingEDA::~RankingEDA()
{
    delete m_best;
    delete m_population;
    delete m_model;
}


/*
 * Running function
 */
int RankingEDA::Run(){
    
    //variable initializations.
    double newScore;
    int * genes= new int[m_problem_size];
    int iterations=1;
    double rate=0.1;
#ifdef PRINT_BEST
    ofstream output_file_log;
    output_file_log.open(m_log_filename);
#endif
#ifdef PRINT_AVG_THETAS
    ofstream output_file_thetas;
    output_file_thetas.open(m_thetas_log_filename);
#endif
    
    //EDA iteration. Stopping criterion is the maximum number of evaluations performed
    while (m_evaluations<m_max_evaluations){

        //learn model
        m_model->Learn(m_population, m_sel_size);
        
        //sample de model
        int accepted=0;
        if ((m_evaluations+m_offspring_size)<m_max_evaluations){
            accepted+=m_model->Sample(m_population,m_offspring_size, m_inverse,m_problem);
            m_evaluations+= m_offspring_size;
        }
        else{
            accepted+=m_model->Sample(m_population,(int)(m_max_evaluations-m_evaluations), m_inverse,m_problem);
            m_evaluations+=(int)(m_max_evaluations-m_evaluations);
        }
        
        
        //update the model.
        m_population->SortPopulation(1);

        //update indicators
        newScore=m_population->m_individuals[0]->Value();
        #if PRINT_AVG_THETAS==1
        double * array=((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_theta_parameters;
        double avg=Average(array,m_problem_size-1);
        output_file_thetas<<avg/GetUpperThetaBound(m_problem_size)<<endl;
        #endif
        
        
#ifdef PRINT_BEST
        output_file_log<<setprecision(10)<<m_population->m_individuals[0]->Value()<<endl;
#endif

        if (newScore>m_best->Value())
        {
            m_best->SetGenes(m_population->m_individuals[0]->Genes());
            m_best->SetValue(newScore);
            m_convergence_evaluations=m_evaluations;
            #if PRINT_LOG==1
           if (((string)m_model_type)=="M")
                cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<"  Theta. "<<((CMallowsModel*)m_model)->m_distance_model->m_theta_parameter<<" lower: "<<m_lower_theta_bound <<" upper: "<<m_upper_theta_bound<<endl;
           else if (((string)m_model_type)=="GM")
                cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<"  Thetas. "<<((CGeneralizedMallowsModel*)m_model)->m_distance_model->m_theta_parameters[0]<<" lower: "<<m_lower_theta_bound <<" upper: "<<m_upper_theta_bound<<endl;
           else {
               cout<<""<<m_population->m_individuals[0]->Value()<<" , "<<m_evaluations<<" , "<<m_max_evaluations-m_evaluations<<endl;
           }

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
    delete [] genes;
#ifdef PRINT_AVG_THETAS
    output_file_thetas.close();
#endif
#ifdef PRINT_BEST
    output_file_log.close();
#endif
    return 0;
}

int RankingEDA::CalculateDistanceAndX(int * sigma, int *x, int m_problem_size){
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
double RankingEDA::CalculateAverageDistace_InPopulation_(CPopulation * population, int sel_size){
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
 * Returns the number of performed evaluations.
 */
long int RankingEDA::GetPerformedEvaluations(){
    return m_convergence_evaluations;
}

/*
 * Returns the fitness of the best solution obtained.
 */
double RankingEDA::GetBestSolutionFitness(){
    return m_best->Value();
}

/*
 * Returns the best solution obtained.
 */
CIndividual * RankingEDA::GetBestSolution(){
    return m_best;
}


/*
 * This method applies a swap of the given i,j positions in the array.
 */
void RankingEDA::Swap(int * array, int i, int j)
{
	int aux=array[i];
	array[i]=array[j];
	array[j]=aux;
}


/*
 * Returns the appropriate upper bound under the KEndall distance for the specified problem size.
 */
double RankingEDA::GetUpperThetaBound_Kendall(int problem_size){
    double upper;
    switch (problem_size) {
        case 10:
            upper=4.4;
            break;
        case 20:
            upper=5.2;
            break;
        case 30:
            upper=5.5;
            break;
        case 40:
            upper=5.8;
            break;
        case 50:
            upper=6.1;
            break;
        case 60:
            upper=6.2;
            break;
        case 70:
            upper=6.3;
            break;
        case 80:
            upper=6.4;
            break;
        case 90:
            upper=6.6;
            break;
        case 100:
            upper=6.7;
            break;
        case 150:
            upper=7.2;
            break;
        case 200:
            upper=7.5;
            break;
        case 250:
            upper=7.7;
            break;
        case 300:
            upper=8.0;
            break;
        case 350:
            upper=8.1;
            break;
        case 400:
            upper=8.2;
            break;
        case 450:
            upper=8.3;
            break;
        case 500:
            upper=8.4;
            break;
        default:
            upper=14;
            break;
    }
    return upper;
}

/*
 * Returns the appropriate upper bound under the Cayley distance for the specified problem size.
 */
double RankingEDA::GetUpperThetaBound_Cayley(int problem_size){
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

