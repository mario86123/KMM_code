#!/bin/bash
# ProblemSizeArray=("20" "40" "60")
ProblemSizeArray=("76")
Problem=("pr76")
# ProblemSizeArray=("60")
PopulationSizeArray=("76" "228" "380" "760" "1520")
SelectionPressureArray=("2" "5" "10" "15" "20")
# SelectionPressureArray=("15" "20")


# ell
for ((ell = 0; ell < 1; ell++))
do

  # Population Size
  for ((pop = 0; pop < 5; pop++))
  do

    # selection pressure
    for ((s = 0; s < 5; s++))
    do

  
      # # repeat times
      for ((times = 0; times < 10; times++))
      do
        
        
        # KMC
        # Problem
        for ((problem = 0; problem < 1; problem++))
        do
        
            ./RankingEDAsCEC -i ./TSP_instance/${Problem[$problem]}.tsp\
                              -o ./TSP_result/KMC/${Problem[$problem]}_"$((${PopulationSizeArray[$pop]}))"_"${SelectionPressureArray[$s]}"_KMC_"$times".txt\
                              -s "$times"\
                              -t TSP\
                              -m M -d C -v 0 -x 2\
                              -r "${SelectionPressureArray[$s]}"\
                              -z "$((${PopulationSizeArray[$pop]}))" &

        done
        # problem end


      done
      # repeat times end

      
    done
    # selection pressure end
    echo "end of pop_size = "$((${PopulationSizeArray[$pop]}))", Emax = "$((${ProblemSizeArray[$ell]} * ${ProblemSizeArray[$ell]} * 1000))", model = KMC"
    wait




  done
   # Population Size end

done
# ell end