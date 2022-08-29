/*
 *  LocalSearch.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 9/15/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "LocalSearch.h"
#include <iomanip>
#include <algorithm>
#include "LocalOptima.h"
#include <sys/time.h>
#include <cmath>
#define MAX(A,B) ( (A > B) ? A : B)

/*
 * The constructor.
 */
LocalSearch::LocalSearch(PBP * problem, int size, int neighborhood, double max_evaluations){
    m_problem=problem;
    m_problem_size=size;
    m_neighborhood=neighborhood;
    m_max_evaluations=max_evaluations;
    m_evaluations=0;
    m_perm_aux= new int[m_problem_size];
	m_best= new int[m_problem_size];
    m_genes_aux=new int[m_problem_size];
}

/*
 * The destructor.
 */
LocalSearch::~LocalSearch(){
    delete m_problem;
    delete [] m_genes_aux;
	delete [] m_perm_aux;
	delete [] m_best;
}

/*
 * Executes a greedy local search for the specified neighborhood until the obtained solution is a local optima.
 */
double LocalSearch::Run(int * solution, LocalOptima * lo){
    
        double fitness,cost_opt,cost_neigh;
        int  * aux_perm=new int[m_problem_size];
        int seguir,j,k,aux;
        int best_k, best_j;
        fitness=m_problem->Evaluate(solution);
        cost_opt=fitness;
        memcpy(aux_perm, solution, sizeof(int)*m_problem_size);
        
        seguir=1;
        switch (m_neighborhood) {
            case 1:
                while (seguir) {
                    seguir=0;
                    for (k=0; k<m_problem_size-1; k++) {
                        aux=aux_perm[k];
                        aux_perm[k]=aux_perm[k+1];
                        aux_perm[k+1]=aux;
                        
                        cost_neigh=m_problem->Evaluate(aux_perm);
                        if (cost_neigh>cost_opt) {

                            best_k=k;
                            cost_opt=cost_neigh;
                            seguir=1;
                        }
                        aux= aux_perm[k+1];
                        aux_perm[k+1]=aux_perm[k];
                        aux_perm[k]=aux;
                    }
                    if (seguir) {
                        aux=aux_perm[best_k];
                        aux_perm[best_k]=aux_perm[best_k+1];
                        aux_perm[best_k+1]=aux;
                    }
                }
                break;
                
            case 2:
                while (seguir) {
                    seguir=0;
                    for (k=0; k<m_problem_size-1; k++) {
                        for (j=k+1; j<m_problem_size; j++) {
                            aux=aux_perm[k];
                            aux_perm[k]=aux_perm[j];
                            aux_perm[j]=aux;
                            
                            cost_neigh=m_problem->Evaluate(aux_perm);
                            if (cost_neigh>cost_opt) {
                                
                                best_k=k;
                                best_j=j;
                                cost_opt=cost_neigh;
                                seguir=1;
                            }
                            aux= aux_perm[j];
                            aux_perm[j]=aux_perm[k];
                            aux_perm[k]=aux;
                        }
                    }
                    if (seguir) {
                        aux=aux_perm[best_k];
                        aux_perm[best_k]=aux_perm[best_j];
                        aux_perm[best_j]=aux;
					}
				}
				break;
                
            case 3:
                int cost_improve, cost_best_provisional;
                int i,j;
                
                memcpy(m_best,aux_perm, sizeof(int)*m_problem_size);
                memcpy(m_genes_aux,aux_perm,sizeof(int)*m_problem_size);
                cost_opt=m_problem->Evaluate(solution);
                cost_best_provisional=cost_opt;
                
                seguir=1;
                while (seguir)
                {
                    seguir=0;
                    //forward step
                    for (i=0;i<(m_problem_size-1);i++)
                    {
                        memcpy(m_perm_aux,m_genes_aux, sizeof(int)*m_problem_size);
                        Swap(m_perm_aux,i,i+1);
                        cost_improve=m_problem->Evaluate(m_perm_aux);
                        if (cost_improve>cost_best_provisional)
                        {
                            memcpy(m_best, m_perm_aux,sizeof(int)*m_problem_size);
                            cost_best_provisional=cost_improve;
                            seguir=1;
                        }
                        
                        for (j=i+1; j<(m_problem_size-1);j++)
                        {
                            Swap(m_perm_aux,j,j+1);
                            cost_improve=m_problem->Evaluate(m_perm_aux);
                            if (cost_improve>cost_best_provisional)
                            {
                                memcpy(m_best, m_perm_aux,sizeof(int)*m_problem_size);
                                cost_best_provisional=cost_improve;
                                seguir=1;
                            }
                        }
                    }
                    
                    //backward step
                    Swap(m_genes_aux,0,1);
                    for (i=2;i<m_problem_size;i++)
                    {
                        Swap(m_genes_aux,0,i);
                        memcpy(m_perm_aux, m_genes_aux, sizeof(int)*m_problem_size);
                        cost_improve=m_problem->Evaluate(m_perm_aux);
                        if (cost_improve>cost_best_provisional)
                        {
                            memcpy(m_best, m_perm_aux,sizeof(int)*m_problem_size);
                            cost_best_provisional=cost_improve;
                            seguir=1;
                        }
                        for (j=0;j<(i-2);j++)
                        {
                            Swap(m_perm_aux,j,j+1);
                            cost_improve=m_problem->Evaluate(m_perm_aux);
                            if (cost_improve>cost_best_provisional)
                            {
                                memcpy(m_best, m_perm_aux,sizeof(int)*m_problem_size);
                                cost_best_provisional=cost_improve;
                                seguir=1;
                            }
                        }
                    }
                    if (seguir==1)
                    {
                        memcpy(m_genes_aux, m_best, sizeof(int)*m_problem_size);
                        cost_opt=cost_best_provisional;
                    }
                }
                memcpy(aux_perm, m_best, sizeof(int)*m_problem_size);
                break;
                
			default:
				break;
        }
    
    memcpy(lo->genes,aux_perm,sizeof(int)*m_problem_size);
    lo->fitness=cost_opt;
    lo->attraction_basin=1;
    
    delete []aux_perm;
    return cost_opt;
}


/*
 * Executes a greedy local search for the specified neighborhood until the obtained solution is a local optima. Besides, the maximum number of evalutions is not exceeded.
 */
double LocalSearch::Run_Flagged(int * solution, LocalOptima * lo){
    
    double cost_opt;
    int  * aux_perm=new int[m_problem_size];
    memcpy(aux_perm, solution, sizeof(int)*m_problem_size);
    
/*    switch (m_neighborhood) {
        case 1:
            cost_opt=LocalSearch_AdjacentSwap(m_problem, aux_perm, m_problem_size);
            break;
            
        case 2:*/
            cost_opt=LocalSearch_Swap(m_problem, aux_perm, m_problem_size);
/*            break;
            
        case 3:
            cost_opt=LocalSearch_Insert(m_problem, aux_perm, m_problem_size);
            break;
            
        default:
            break;
    }
  */  
    memcpy(lo->genes,aux_perm,sizeof(int)*m_problem_size);
    lo->fitness=cost_opt;
    lo->attraction_basin=1;
    
    delete []aux_perm;
    return cost_opt;
}

/*
 * This method applies a greedy local search with the insert operator neighborhood.
 */
double LocalSearch::LocalSearch_Insert(PBP * problem, int * solution, int size){
    
	int cost_opt, cost_improve, cost_best_provisional;
	int seguir, i,j;

    memcpy(m_best,solution, sizeof(int)*size);
    memcpy(m_genes_aux,solution,sizeof(int)*size);
    m_evaluations++;
	cost_opt=problem->Evaluate(solution);
	cost_best_provisional=cost_opt;
	//double variation=0;
    //double ratio=0;
    
	seguir=1;
	//int iter=0;
	while (seguir && m_evaluations<m_max_evaluations)
	{
		seguir=0;
		//forward step
		for (i=0;i<(size-1) && m_evaluations<m_max_evaluations;i++)
		{
            memcpy(m_perm_aux,m_genes_aux, sizeof(int)*size);
			Swap(m_perm_aux,i,i+1);
			cost_improve=problem->Evaluate(m_perm_aux);
            m_evaluations++;
          //  variation+=abs(cost_improve-cost_best_provisional);
			if (cost_improve>cost_best_provisional)
			{
                memcpy(m_best, m_perm_aux,sizeof(int)*size);
				cost_best_provisional=cost_improve;
				seguir=1;
			}
			
			for (j=i+1; j<(size-1) && m_evaluations<m_max_evaluations;j++)
			{
				Swap(m_perm_aux,j,j+1);
				cost_improve=problem->Evaluate(m_perm_aux);
                m_evaluations++;
          //  variation+=abs(cost_improve-cost_best_provisional);
				if (cost_improve>cost_best_provisional)
				{
                    memcpy(m_best, m_perm_aux,sizeof(int)*size);
					cost_best_provisional=cost_improve;
					seguir=1;
				}
			}
		}
        
		//backward step
		Swap(m_genes_aux,0,1);
		for (i=2;i<size && m_evaluations<m_max_evaluations;i++)
		{
			Swap(m_genes_aux,0,i);
			memcpy(m_perm_aux, m_genes_aux, sizeof(int)*size);
            cost_improve=problem->Evaluate(m_perm_aux);
            m_evaluations++;
            //variation+=abs(cost_improve-cost_best_provisional);
			if (cost_improve>cost_best_provisional)
			{
                memcpy(m_best, m_perm_aux,sizeof(int)*size);
				cost_best_provisional=cost_improve;
				seguir=1;
			}
			for (j=0;j<(i-2) && m_evaluations<m_max_evaluations;j++)
			{
				Swap(m_perm_aux,j,j+1);
				cost_improve=problem->Evaluate(m_perm_aux);
                m_evaluations++;
               // variation+=abs(cost_improve-cost_best_provisional);
				if (cost_improve>cost_best_provisional)
				{
                    memcpy(m_best, m_perm_aux,sizeof(int)*size);
					cost_best_provisional=cost_improve;
					seguir=1;
				}
			}
		}
        //ratio+=variation/(m_problem_size*m_problem_size - 2*m_problem_size +1);
      //  variation=0;
		if (seguir==1)
		{
            memcpy(m_genes_aux, m_best, sizeof(int)*size);
			cost_opt=cost_best_provisional;
		}
	//	iter++;
	}
    //cout<<setprecision(10)<<ratio/iter<<endl;
    memcpy(solution, m_best, sizeof(int)*size);

	return cost_opt;
}

/*
 * Calculates the attraction basin.
 */
vector<LocalOptima*> LocalSearch::Calculate_AttBasin(){
    
    vector<LocalOptima*> localoptima;
    int numlocaloptima=0;

    int iter=0;
    int j;
    int * solution= new int[m_problem_size];
    for (int i=0;i<m_problem_size;i++) solution[i]=i;
    sort (solution,solution+m_problem_size);
    LocalOptima * lo= new LocalOptima(m_problem_size);
    do{
       // PrintArray(solution, m_problem_size, "");
        (this)->Run(solution,lo);
        
        for (j=0;j<localoptima.size();j++)
        {
            if (localoptima[j]->Equal(lo, m_problem_size)){
                //incremente attbasin size
                localoptima[j]->attraction_basin++;
                break;
            }
        }
        if (j==localoptima.size()){
            //add at the end new local optima
            numlocaloptima++;
            localoptima.resize(numlocaloptima);
            localoptima[numlocaloptima-1]=new LocalOptima(m_problem_size,lo);
        }

        iter++;
    }
    while ( next_permutation (solution,solution+m_problem_size) );
    sort(localoptima.begin(), localoptima.end(), Better);
    delete [] solution;
    
    return localoptima;
}

/*
 * This method applies a greedy local search with the swap operator neighborhood.
 */
double LocalSearch::LocalSearch_Swap(PBP * problem, int * solution, int size){
    double cost_opt,cost_neigh;
    int seguir,j,k,aux;
    int best_k, best_j;
    cost_opt=m_problem->Evaluate(solution);
    m_evaluations++;

    while (seguir && m_evaluations<m_max_evaluations) {
        seguir=0;
        for (k=0; k<m_problem_size-1 && m_evaluations<m_max_evaluations; k++) {
            for (j=k+1; j<m_problem_size && m_evaluations<m_max_evaluations; j++) {
                aux=solution[k];
                solution[k]=solution[j];
                solution[j]=aux;
                
                cost_neigh=m_problem->Evaluate(solution);
                m_evaluations++;
                //      variation+=abs(cost_neigh-cost_opt);
                if (cost_neigh>cost_opt) {
                    
                    best_k=k;
                    best_j=j;
                    cost_opt=cost_neigh;
                    seguir=1;
                }
                aux= solution[j];
                solution[j]=solution[k];
                solution[k]=aux;
            }
        }
        //ratio+=variation/(m_problem_size*(m_problem_size-1)/2);
        //   iterations++;
        if (seguir) {
            //      variation=0;
            aux=solution[best_k];
            solution[best_k]=solution[best_j];
            solution[best_j]=aux;
            //  cost_opt=m_problem->Evaluate(aux_perm);
        }
    }
    // cout<<setprecision(15)<<ratio/iterations<<endl;*/
	return cost_opt;
}

/*
 * This method applies a greedy local search with the adjacent swap operator neighborhood.
 */
double LocalSearch::LocalSearch_AdjacentSwap(PBP * problem, int * solution, int size){
    int seguir=1;
    double cost_opt,cost_neigh;
    int k,aux;
    int best_k;
    cost_opt=m_problem->Evaluate(solution);
    m_evaluations++;
    
    while (seguir && m_evaluations<m_max_evaluations) {
        seguir=0;
        for (k=0; k<m_problem_size-1 && m_evaluations<m_max_evaluations; k++) {
            aux=solution[k];
            solution[k]=solution[k+1];
            solution[k+1]=aux;
            
            cost_neigh=m_problem->Evaluate(solution);
            m_evaluations++;
            // variation+=abs(cost_neigh-cost_opt);
            if (cost_neigh>cost_opt) {
                
                best_k=k;
                cost_opt=cost_neigh;
                seguir=1;
            }
            aux= solution[k+1];
            solution[k+1]=solution[k];
            solution[k]=aux;
        }
        //   ratio+=variation/(m_problem_size-1);
        // iterations++;
        if (seguir) {
            //       variation=0;
            aux=solution[best_k];
            solution[best_k]=solution[best_k+1];
            solution[best_k+1]=aux;
            //  cost_opt=m_problem->Evaluate(aux_perm);
        }
    }
    //     cout<<setprecision(15)<<ratio/iterations<<endl;
  
	return cost_opt;
}



