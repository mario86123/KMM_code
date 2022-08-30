#!/bin/bash
# ProblemSizeArray=("20" "40" "60")
ProblemSizeArray=("100")
Problem=("100_5" "100_10" "100_20")
# ProblemSizeArray=("60")
PopulationSizeArray=("3" "5" "10")
SelectionPressureArray=("2" "5" "10")


# ell
for ((ell = 0; ell < 1; ell++))
do

  # Population Size
  for ((pop = 0; pop < 3; pop++))
  do

    # selection pressure
    for ((s = 0; s < 3; s++))
    do

  
      # # repeat times
      for ((times = 0; times < 10; times++))
      do
        
        
        # KMC
        # Problem
        for ((problem = 0; problem < 3; problem++))
        do
        
            ./RankingEDAsCEC -i ./PFSP_instance/tai${Problem[$problem]}.txt\
                              -o ./PFSP_result/KMC/${Problem[$problem]}_"$((${PopulationSizeArray[$pop]} * ${ProblemSizeArray[$ell]}))"_"${SelectionPressureArray[$s]}"_KMC_"$times".txt\
                              -s "$times"\
                              -t PFSP\
                              -m M -d C -v 0 -x 2\
                              -r "${SelectionPressureArray[$s]}"
                              -x "$((${PopulationSizeArray[$pop]} * ${ProblemSizeArray[$ell]}))" &

            echo "end of problem = tai${Problem[$problem]}.txt, pop_size = "$((${PopulationSizeArray[$pop]} * ${ProblemSizeArray[$ell]}))", Emax = "$((${ProblemSizeArray[$ell]} * ${ProblemSizeArray[$ell]} * 1000))", model = GMC"

        done
        # problem end


      done
      # repeat times end
      wait

      
    done
    # selection pressure end



  done
   # Population Size end

done
# ell end