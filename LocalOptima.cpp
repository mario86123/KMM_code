//
//  LocalOptima.cpp
//  MultistartLocalSearch
//
//  Created by Josu Ceberio Uribe on 04/03/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//
#include "LocalOptima.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include "string.h"
using std::cerr;
using std::cout;
using std::endl;

/*
 * The constructor.
 */
LocalOptima::LocalOptima(int size){
    genes= new int[size];
    attraction_basin=0;
    fitness=0;
    problem_size=size;
}

/*
 * The constructor.
 */
LocalOptima::LocalOptima(int size, LocalOptima * lo){
    genes= new int[size];
    problem_size=size;
    attraction_basin=lo->attraction_basin;
    fitness=lo->fitness;
    memcpy(genes,lo->genes, sizeof(int)*size);
}

/*
 * The descructor of the class individual
 */
LocalOptima::~LocalOptima()
{
	delete [] genes;
	genes=NULL;
}

/*
 * Is equal to the local optima.
 */
bool LocalOptima::Equal(LocalOptima * lo, int size){
    return (lo->fitness==fitness && memcmp(lo->genes,genes,sizeof(int)*size)==0);
}

/*
 * Prints the local optima solution.
 */
void LocalOptima::Print(){
    cout<<" "<<fitness<<"     "<<attraction_basin<<"        [";
    for (int i=0;i<problem_size;i++){
        cout<<" "<<genes[i];
    }
    cout<<" ]"<<endl;
}

/*
 * Prints the local optima solution.
 */
void LocalOptima::PrintFile(ofstream file){
    file<<" "<<fitness<<"; "<<attraction_basin<<"; [";
    for (int i=0;i<problem_size;i++){
        file<<" "<<genes[i];
    }
    file<<" ]"<<endl;
}

ostream & operator<<(ostream & os,LocalOptima * & lo)
{
    os<< lo->fitness<<"; "<<lo->attraction_basin<<"; [";
	for(int i=0;i<lo->problem_size;i++)
        os << " " << lo->genes[i];
    os<<" ]";
	return os;
}

