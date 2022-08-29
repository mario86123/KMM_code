/*
 *  Tools.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <vector>

using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;
using std::string;

/*
 * Returs the first position at which value appears in array. If it does not appear, then it returns -1;
 */
int Find(int * array, int size, int value);

/*
 * Calculates the average value of the parameters in the array.
 */
double Average(double * params, int size);

double Variance(double * params, int size, double avg );

/*
 * Calculates Kullback-Leibeler divergence between P and Q distributions.
 */
double KullbackLeibelerDivergence(double * P, double * Q, int size);

/*
 * Calculates Total Variation divergence between P and Q distributions.
 */
double TotalVariationDivergence(double * P, double * Q, int size);

/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int * permutation, int size);

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int * permutation, int n);

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr);

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int* array, int length, string text);

/*
 * Prints in standard output 'length' long double elements of a given array.
 */
void PrintArray(long double* array, int length, string text);

/*
 * Prints the given doubles array in the standard output.
 */
void PrintArray(double* array, int length, string text);

/*
 * Prints the given integers matrix in the standard output.
 */
void PrintMatrix(int** matrix, int length, int length2, string text);

/*
 * Prints the given doubles matrix in the standard output.
 */
void PrintMatrix(double** matrix, int length, int length2, string text);

/*
 * Applies the random keys sorting strategy to the vector of doubles
 */
void RandomKeys( int * a, double * criteriaValues, int size);

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size);

/*
 * Calculates the Kendall tau distance between 2 permutations.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux);

/*
 * Calculates the Kendall tau distance between 2 permutations.
 * Auxiliary parameters are used for multiple executions.
 */
int Kendall(int* permutationA, int*permutationB, int size, int * m_aux, int * invertedB, int *  composition, int * v);

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size);

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int * permutationA, int * permutationB, int size,int * invertedB, int *  composition, int * elemsToCycles, int * maxPosInCycle, int * freeCycle);
int FindNewCycle(int * freeCycle, int size);
int NextUnasignedElem(int * elemsToCycles, int size);
int CalculateDistance(int*sigma, int size);


/*
 * Calculates the length of the longest increasing subsequence in the given array of ints.
 */
int getLISLength(int*sigma,int size);

/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int*s1, int*s2, int*res, int n);

/*
* Calculates V_j-s vector.
*/
void vVector(int*v, int*permutation, int n);

/*
 *  Optimized version by Leti of the V_j-s vector calculation.
 */
void vVector_Fast(int*v, int*permutation, int n, int * m_aux);

/*
 * Inverts a permutation.
 */
void Invert(int*permu, int n, int* inverted);

/*
 * This method moves the value in position i to the position j.
 */
void InsertAt(int * array, int i, int j, int n);

/*
 * Calculates the factorial of a solution.
 */
long double factorial(int val);

/*
 * This method applies a swap of the given i,j positions in the array.
 */
void Swap(int * array, int i, int j);


