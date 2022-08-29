/*
 *  Tools.cpp
 *  DiscreteEDA
 *
 *  Created by Josu Ceberio Uribe on 9/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */
#include "PermutationTools.h"
#include "Tools.h"
#include <string.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using std::stringstream;
using std::string;
using std::cerr;
using std::cout;
using std::endl;

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int n, int * permutation) 
{
	for (int i = 0; i < n; ++i) 
	{
		int j = rand() % (i + 1);
		permutation[i] = permutation[j];
		permutation[j] = i;
	}
}
