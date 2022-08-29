//
//  LocalOptima.h
//  MultistartLocalSearch
//
//  Created by Josu Ceberio Uribe on 04/03/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __MultistartLocalSearch__LocalOptima__
#define __MultistartLocalSearch__LocalOptima__
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "LocalOptima.h"
using std::istream;
using std::ostream;
using std::string;
using namespace std;
class LocalOptima
{
public:
    

    /*
     * The genes of the solution.
     */
    int * genes;
    
    /*
     * The fitness of the solution.
     */
    double fitness;
    
    /*
     * The size of the attraction basin.
     */
    double attraction_basin;
    
    /*
     * The size of the problem.
     */
    int problem_size;
    /*
     * The constructor.
     */
	LocalOptima(int length);
    
    /*
     * The constructor.
     */
    LocalOptima(int size, LocalOptima * lo);
    
	/*
     * The destructor.
     */
	virtual ~LocalOptima();
    
    /*
     * Is equal to the local optima.
     */
    bool Equal(LocalOptima * lo, int size);
    
    /*
     * Prints the local optima solution.
     */
    void Print();
    
    /*
     * Prints the local optima solution.
     */
    void PrintFile(ofstream results_file);
};

#endif /* defined(__MultistartLocalSearch__LocalOptima__) */
