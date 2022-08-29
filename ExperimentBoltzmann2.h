//
//  ExperimentBoltzmann2.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 10/11/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__ExperimentBoltzmann2__
#define __RankingEDAsCEC__ExperimentBoltzmann2__

#include <stdio.h>
#include "LOP.h"
class CExperimentBoltzmann2
{
    
public:
    
    /*
     * The constructor.
     */
    CExperimentBoltzmann2();
    
    /*
     * The destructor.
     */
    virtual ~CExperimentBoltzmann2();
    
    /*
     * Running function.
     */
    void Run(int problem_size, LOP * problem);

private:
    
};

#endif /* defined(__RankingEDAsCEC__ExperimentBoltzmann2__) */
