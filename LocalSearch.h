/*
 *  LocalSearch.h
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 9/15/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */
#ifndef _LOCALSEARCH_H__
#define _LOCALSEARCH_H__

#include <fstream>
#include <iostream>
#include "LocalOptima.h"
#include "Tools.h"
#include "PBP.h"
using std::istream;
using std::ostream;
using std::string;
using namespace std;
class LocalSearch
{
	
public:

    struct Better{
        bool operator()(LocalOptima* a, LocalOptima * b) const{
            return a->fitness>b->fitness;
        }
    } Better;
    
    int * m_perm_aux;
	int * m_best;
    int * m_genes_aux;
    /*
     * Problem
     */
    PBP * m_problem;
    
    /*
     * Problem size
     */
    int m_problem_size;
    
    /*
     * Maximum number of evaluations permitted in the execution.
     */
    double m_max_evaluations;
    
    /*
     * The number of evaluations performed at the moment.
     */
    double m_evaluations;
    
    /*
     * The neighborhood.
     */
    int m_neighborhood;
    
    /*
     * The constructor.
     */
    LocalSearch(PBP * problem, int size, int neighborhood, double max_evaluations);
    
    /*
     * The destructor.
     */
    virtual ~LocalSearch();
    
    
    /*
     * Executes a greedy local search for the specified neighborhood.
     */
    double Run(int * solution, LocalOptima * lo);
    
    /*
     * Executes a greedy local search for the specified neighborhood until the obtained solution is a local optima. Besides, the maximum number of evalutions is not exceeded.
     */
    double Run_Flagged(int * solution, LocalOptima * lo);
    
    /*
     * Calculates the attraction basin.
     */
    vector<LocalOptima*> Calculate_AttBasin();

private:
    
    /*
     * This method applies a greedy local search with the insert operator neighborhood.
     */
    double LocalSearch_Insert(PBP * problem, int * solution, int size);
    
    /*
     * This method applies a greedy local search with the swap operator neighborhood.
     */
    double LocalSearch_Swap(PBP * problem, int * solution, int size);
    
    /*
     * This method applies a greedy local search with the adjacent swap operator neighborhood.
     */
    double LocalSearch_AdjacentSwap(PBP * problem, int * solution, int size);

    
};
#endif
