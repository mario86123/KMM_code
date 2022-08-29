//
//  ExperimentBoltzmann.h
//  RankingEDAsCEC
//
//  Created by Josu Ceberio Uribe on 04/11/14.
//  Copyright (c) 2014 Josu Ceberio Uribe. All rights reserved.
//

#ifndef __RankingEDAsCEC__ExperimentBoltzmann__
#define __RankingEDAsCEC__ExperimentBoltzmann__

#include <stdio.h>
#include "PBP.h"

class CExperimentBoltzmann
{
    
public:
 
    /*
     * The constructor.
     */
    CExperimentBoltzmann();
    
    /*
     * The destructor.
     */
    virtual ~CExperimentBoltzmann();
    
    /*
     * Running function.
     */
    void Run(int problem_size, PBP * problem);
    
private:
    
};
#endif /* defined(__RankingEDAsCEC__ExperimentBoltzmann__) */
