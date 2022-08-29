/*
 *  LOP.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "LOP.h"
#include "Tools.h"

/*
 * Class constructor.
 */
LOP::LOP()
{
	
    
}

/*
 * Class destructor.
 */
LOP::~LOP()
{
	for (int i=0;i<m_problemsize;i++)
		delete [] m_matrix[i];
	delete [] m_matrix;
    delete [] m_aux;
}

/*
 * Read LOP instance file.
 */
int LOP::Read(string filename)
{
	char line[2048]; // variable for input value
	string data="";
	ifstream indata;
	indata.open(filename.c_str(),ios::in);
	int num=0;
	while (!indata.eof())
	{

		indata.getline(line, 2048);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
		if (sline=="")
		{
			break;
		}
		if (num==0)
		{
			m_problemsize = atoi(line);
		}
		else
		{
			if (data=="")
				data = line;
			else
				data = data+' '+line;
		}
		num++;
	}
	indata.close();

	//BUILD MATRIX
	m_matrix = new double*[m_problemsize];
	for (int i=0;i<m_problemsize;i++)
	{
		m_matrix[i]= new double[m_problemsize];
	}
    m_aux= new int[m_problemsize];
    
	istringstream iss(data);
	int i=0;
	int j=0;
	do
	{
		string sub;
	    iss >> sub;
	    if (sub!=""){
			//save distance in distances matrix.
	    	m_matrix[i][j]= atof(sub.c_str());
	    	if (j==(m_problemsize-1))
	    	{
	    		i++;
	    		j=0;
	    	}
	    	else
	    	{
	    		j++;
	    	}
	    }
	    else
	    {
	    	break;
	    }
	} while (iss);

	return (m_problemsize);
}

/*
 * This function evaluates the inverted solution of the given solution for the LOP problem.
 */
double LOP::EvaluateInv(int * genes)
{
    Invert(genes,m_problemsize,m_aux); //original
    double fitness=0;

	
    int i,j;
	for (i=0;i<m_problemsize-1;i++)
		for (j=i+1;j<m_problemsize;j++)
			fitness+= m_matrix[m_aux[i]][m_aux[j]];
	return fitness;
}

/*
 * This function evaluates the solution for the LOP problem.
 */
double LOP::Evaluate(int * genes)
{
	double fitness=0;
    int i,j;
	for (i=0;i<m_problemsize-1;i++)
		for (j=i+1;j<m_problemsize;j++)
			fitness+=(double) m_matrix[genes[i]][genes[j]];
	return fitness;
}

/*
 * Returns the size of the problem.
 */
int LOP::GetProblemSize()
{
    return m_problemsize;
}

/*
 * Calculates the exact Boltzmann distribution for the LOP.
 */
double * LOP::Calculate_BoltzmannDistribution(double c){
    
  /*  int num_permus=factorial(m_problemsize);
    double * boltzmann_distribution= new double[num_permus];

    //calculate the normalization constant by bruteforce.
    double Z_2=0;
    int p1[m_problemsize];
    for (int i=0;i<m_problemsize;i++){
        p1[i]=i;
    }
    sort(p1,p1+m_problemsize);

    double eval=0;
    int index=0;
    do{
        eval=(double)Evaluate(p1)/c;
        Z_2 += exp(eval);
        boltzmann_distribution[index] =  exp(eval);
        index++;
    }
    while ( next_permutation (p1,p1+m_problemsize) );

    //calculate the boltzmann distribution.
    for (int i=0;i<num_permus;i++){
        boltzmann_distribution[i]=boltzmann_distribution[i]/Z_2;
        
    }
    return boltzmann_distribution;*/
    return 0;
}

double LOP::Contribution(int * solution, int item, int position){
    double contribution=0;
    int i;
    for (i=0;i<position;i++){
        contribution+=m_matrix[solution[i]][item];

    }
    for (i=position+1;i<m_problemsize;i++){
        contribution+=m_matrix[item][solution[i]];
    }
//    cout<<contribution<<"item: "<<item<<" pos: "<<position<<endl;

    return contribution;
}









