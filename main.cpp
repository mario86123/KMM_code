//
//  main.cpp
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 11/19/13.
//  Copyright (c) 2013 Josu Ceberio Uribe. All rights reserved.
//
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include "RankingModel.h"
#include "MallowsModel.h"
#include "PlackettLuce.h"
#include "RankingEDA.h"
#include "ExperimentBoltzmann2.h"
#include "Distance.h"
#include "MixtureRankingEDA.h"
#include "KernelRankingEDA.h"
#include "Ulam2.h"
#include "PBP.h"
#include "PFSP.h"
#include "LOP.h"
#include "QAP.h"
#include "TSP.h"
#include "API.h"


//It is the instance of the problem to optimize.
PBP * PROBLEM;

//The type of the problem to solve.
char PROBLEM_TYPE[10];

//The type of the model to use in the EDA.
char MODEL_TYPE[10];

//The type of the metric to use in the model of the EDA.
char METRIC_TYPE[10];

// The individual size.
int PROBLEM_SIZE;

// Number of evaluations performed.
long int EVALUATIONS = 0;

// Convergence evaluation of the best fitness.
long int CONVERGENCE_EVALUATIONS = 0;

// Maximum number of evaluations allowed performed.
long int MAX_EVALUATIONS = 0;

// Name of the file where the result will be stored.
char RESULTS_FILENAME[50];

// Name of the file where the instances is stored.
char INSTANCE_FILENAME[50];

// Name of the file where the logs will be stored.
char LOG_FILENAME[50];

// Name of the file where the logs of the thetas will be stored.
char THETAS_LOG_FILENAME[50];

// Name of the file where the logs of the weights will be stored.
char WEIGHTS_LOG_FILENAME[50];

// The seed asigned to the process
int SEED;

//Determines if a solution needs to be inverted before evaluating.
int INVERSE;

//Determines if mixtures of models are learnt, kernels or on the contrary, only one model is learnt.
int EDA_TYPE;

//The number of clusters to learn in the mixture.
int NUM_CLUSTERS;


/*
 * Get next command line option and parameter
 */
int GetOption (int argc, char** argv, char* pszValidOpts, char** ppszParam)
{
	
    static int iArg = 1;
    char chOpt;
    char* psz = NULL;
    char* pszParam = NULL;
	
    if (iArg < argc)
    {
        psz = &(argv[iArg][0]);
		
        if (*psz == '-' || *psz == '/')
        {
            // we have an option specifier
            chOpt = argv[iArg][1];
			
            if (isalnum(chOpt) || ispunct(chOpt))
            {
                // we have an option character
                psz = strchr(pszValidOpts, chOpt);
				
                if (psz != NULL)
                {
                    // option is valid, we want to return chOpt
                    if (psz[1] == ':')
                    {
                        // option can have a parameter
                        psz = &(argv[iArg][2]);
                        if (*psz == '\0')
                        {
                            // must look at next argv for param
                            if (iArg+1 < argc)
                            {
                                psz = &(argv[iArg+1][0]);
                                if (*psz == '-' || *psz == '/')
                                {
                                    // next argv is a new option, so param
                                    // not given for current option
                                }
                                else
                                {
                                    // next argv is the param
                                    iArg++;
                                    pszParam = psz;
                                }
                            }
                            else
                            {
                                // reached end of args looking for param
                            }
							
                        }
                        else
                        {
                            // param is attached to option
                            pszParam = psz;
                        }
                    }
                    else
                    {
                        // option is alone, has no parameter
                    }
                }
                else
                {
                    // option specified is not in list of valid options
                    chOpt = -1;
                    pszParam = &(argv[iArg][0]);
                }
            }
            else
            {
                // though option specifier was given, option character
                // is not alpha or was was not specified
                chOpt = -1;
                pszParam = &(argv[iArg][0]);
            }
        }
        else
        {
            // standalone arg given with no option specifier
            chOpt = 1;
            pszParam = &(argv[iArg][0]);
        }
    }
    else
    {
        // end of argument list
        chOpt = 0;
    }
	
    iArg++;
	
    *ppszParam = pszParam;
    return (chOpt);
}

/*
 * Help command output.
 */
void usage(char *progname)
{
    cout << "Algorithm -i <instance_name> -o <results_name> -s <seed> -t <problem_type> -m <model_type> -d <metric> -v <inverse> -x <mixture> -k <num_clusters> -l <log_filename> -p <thetas_log_filename> -w <weights_log_filename>" <<endl;
    cout <<"   -i File name of the instance.\n"<<endl;
    cout <<"   -o Name of the file to store the results.\n"<<endl;
    cout <<"   -s Seed to be used for pseudo-random numbers generator.\n"<<endl;
    cout <<"   -t problem_type (TSP, QAP, LOP, PFSP or API).\n"<<endl;
    cout <<"   -m model_type (M (Mallows), GM (Generalized Mallows)).\n"<<endl;
    cout <<"   -d metric (K (Kendall), C (Cayley), U (Ulam)).\n"<<endl;
    cout <<"   -v inverse (0<- no inverse, 1<- inverse)."<<endl;
    cout <<"   -x mixture (0<- normal EDA, 1<- mixtures, 2<- kernels)."<<endl;
    cout <<"   -k number of clusters."<<endl;
    cout <<"   -l Name of the file to store the logs.\n"<<endl;
    cout <<"   -p Name of the file to store the logs of the theta averages.\n"<<endl;
    cout <<"   -w Name of the file to store the logs of the weights of the clusters.\n"<<endl;
}

/*
 * Obtaint the execution parameters from the command line.
 */
bool GetParameters(int argc,char * argv[])
{
	char c;
    if(argc==1)
    {
    	usage(argv[0]);
        return false;
    }
	char** optarg;
	optarg = new char*[argc];
    while ((c = GetOption (argc, argv, "s:h:t:o:i:m:d:v:x:k:l:w:p:",optarg)) != '\0')
    {
    	switch (c)
    	{
                
            case 'h' :
                usage(argv[0]);
                return false;
                break;
                
            case 's' :
                SEED = atoi(*optarg);
                break;
				
           	case 't':
                strcpy(PROBLEM_TYPE, *optarg);
                break;
                
            case 'o' :
                strcpy(RESULTS_FILENAME, *optarg);
                break;
                
			case 'i':
                strcpy(INSTANCE_FILENAME, *optarg);
                break;

            case 'm':
                strcpy(MODEL_TYPE, *optarg);
                break;
                
           	case 'd':
                strcpy(METRIC_TYPE, *optarg);
                break;
                
           	case 'v':
                INVERSE=atoi(*optarg);
                break;
		
            case 'x':
                EDA_TYPE=atoi(*optarg);
                break;
                
            case 'k':
                NUM_CLUSTERS=atoi(*optarg);
                break;

            case 'l':
                strcpy(LOG_FILENAME, *optarg);
                break;
                
            case 'p':
                strcpy(THETAS_LOG_FILENAME, *optarg);
                break;
                
            case 'w':
                strcpy(WEIGHTS_LOG_FILENAME, *optarg);
                break;
        }
    }

    delete [] optarg;
    
	return true;
}

/*
 * Writes the results of the execution.
 */
void WriteResults(double best_fitness, int * best, long int evaluations, long double time_interval){

    ofstream output_file;
    output_file.open(RESULTS_FILENAME);
    output_file<<"Best fitness: "<<setprecision(15)<<best_fitness<<endl;
    output_file<<"Best solution: ";
    for (int i=0;i<PROBLEM_SIZE;i++)
        output_file<<best[i]<<" ";
    output_file<<endl;
    output_file<<"Evaluations performed: "<<setprecision(15)<<evaluations<<endl;
    output_file<<"Time consumed: "<<time_interval<<endl;
    output_file.close();
}

/*
 * Writes the results of the execution.
 */
void WriteResults(double best_fitness, int * best, long int evaluations, long double time_interval, int distance_to_optimum){
    
    ofstream output_file;
    output_file.open(RESULTS_FILENAME);
    output_file<<"Best fitness: "<<setprecision(15)<<best_fitness<<endl;
    output_file<<"Best solution: ";
    for (int i=0;i<PROBLEM_SIZE;i++)
        output_file<<best[i]<<" ";
    output_file<<endl;
    output_file<<"Distance optimum: "<<distance_to_optimum<<endl;
    output_file<<"Evaluations performed: "<<setprecision(15)<<evaluations<<endl;
    output_file<<"Time consumed: "<<time_interval<<endl;
    output_file.close();
}

/*
 * Reads the problem info of the instance set.
 */
PBP * GetProblemInfo(string problemType, string filename)
{
    PBP * problem;
    if (problemType=="PFSP")
        problem= new PFSP();
    else if (problemType=="TSP")
        problem= new TSP();
    else if (problemType=="QAP")
        problem= new QAP();
    else if (problemType=="LOP")
        problem= new LOP();
    else if (problemType=="API")
        problem= new API();
    else{
        cout<<"Wrong problem type was specified."<<endl;
        exit(1);
    }
    
    //Read the instance.
    int problem_size= problem->Read(filename);
    PROBLEM_SIZE=problem_size;
    
	return problem;
}

/*
 * Determines if the problem type of the instance to optimize is an Artifical Permutation Instance (API).
 */
bool IsAPI(string problemType){
    return (problemType=="API");
}

/*
 * Main function.
 */
int main(int argc, char * argv[])
{
   
    //Initialize time variables.
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
    
    //Get parameters
	if(!GetParameters(argc,argv)) return -1;

    //Set seed
    srand(SEED*1000);
    
    //Read the problem instance to optimize.
	PROBLEM = GetProblemInfo(PROBLEM_TYPE,INSTANCE_FILENAME);
    MAX_EVALUATIONS=PROBLEM_SIZE*PROBLEM_SIZE*1000;
        /****************************************************************/
    //Código para sacar los probability charts.
    /*double low=0.0;
    double up=20;
    PROBLEM_SIZE=500;
    cout<<"Problem size: "<<PROBLEM_SIZE<<endl;
    CMallowsModel * model= new CMallowsModel(PROBLEM_SIZE,1,"K",&low,&up);
    int * cr= new int[PROBLEM_SIZE];
    for (int i=0;i<PROBLEM_SIZE;i++) cr[i]=i;
    for (double theta=0.01;theta<20;theta+=0.1){
        model->Learn(cr, theta);
        cout<<model->m_distance_model->Probability(cr)<<endl;
    }
    delete [] cr;
    delete model;
    exit(1);*/
    /****************************************************************/
    //Muestra para el trabajo con Aritz.
   /* int total_num=1000000;
    int print_num=100000;
    CPopulation * pop= new CPopulation(total_num, 1, PROBLEM_SIZE);
    int * perm= new int[PROBLEM_SIZE];
    double fitness;
    for (int i=0;i<total_num;i++)
    {
        GenerateRandomPermutation(perm, PROBLEM_SIZE);
        fitness= PROBLEM->Evaluate(perm);
        pop->SetToPopulation(perm, i,fitness);
    }
    pop->SortPopulation(1);

    for (int i=0;i<print_num;i++){
        PrintArray(pop->m_individuals[i]->Genes(),PROBLEM_SIZE,"");
    }
    delete [] perm;
    delete pop;
    exit(1);*/
    /****************************************************************/
        //Initialize the algorithm
    if (BOLTZMANN_ANALYSIS==1){
        CExperimentBoltzmann2 * experiment= new CExperimentBoltzmann2();
        experiment->Run(PROBLEM_SIZE, (LOP*)PROBLEM);

        delete experiment;
        exit(1);
    }
    long int evaluations;
    CIndividual * best;
    if (EDA_TYPE==1){
        cout<<"Mixture Ranking EDA..."<<endl;
        MixtureRankingEDA * alg= new MixtureRankingEDA(PROBLEM, PROBLEM_SIZE,MAX_EVALUATIONS,MODEL_TYPE,METRIC_TYPE,INVERSE, NUM_CLUSTERS, LOG_FILENAME, THETAS_LOG_FILENAME, WEIGHTS_LOG_FILENAME);
        alg->Run();
        //Get best solution fitness and the number of evaluations performed.
        evaluations = alg->GetPerformedEvaluations();
        best = alg->GetBestSolution();
        
    }
    else if (EDA_TYPE==2){
        cout<<"Kernels Ranking EDA..."<<endl;
        KernelRankingEDA * alg= new KernelRankingEDA(PROBLEM, PROBLEM_SIZE,MAX_EVALUATIONS,MODEL_TYPE,METRIC_TYPE,INVERSE,  LOG_FILENAME, THETAS_LOG_FILENAME);
        alg->Run();
        //Get best solution fitness and the number of evaluations performed.
        evaluations = alg->GetPerformedEvaluations();
        best = alg->GetBestSolution();
        
    }
    else{
        cout<<"Ranking EDA..."<<endl;
        RankingEDA * alg= new RankingEDA(PROBLEM,PROBLEM_SIZE, MAX_EVALUATIONS, MODEL_TYPE, METRIC_TYPE,INVERSE, LOG_FILENAME,THETAS_LOG_FILENAME);
        alg->Run();

        //Get best solution fitness and the number of evaluations performed.
        evaluations = alg->GetPerformedEvaluations();
        best = alg->GetBestSolution();
     
    }
       int dist=0;
        if (IsAPI(PROBLEM_TYPE)==1)
            dist=((API*)PROBLEM)->CalculateDistanceToOptimum(best->Genes());

        gettimeofday(&tim, NULL);
        double t2=tim.tv_sec+(tim.tv_usec/1000000.0);

        //Create the file to store the results.
        ofstream output_file;
        output_file.open(RESULTS_FILENAME);
        output_file.close();
        
        //Print results
        if (IsAPI(PROBLEM_TYPE)==1)
            WriteResults(best->Value(),best->Genes(),evaluations,t2-t1,dist);
        else
            WriteResults(best->Value(),best->Genes(),evaluations,t2-t1);
    
        //delete alg;
   
    return 0;
}

