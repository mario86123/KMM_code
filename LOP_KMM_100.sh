#!/bin/bash
# ProblemSizeArray=("20" "40" "60")
ProblemSizeArray=("100")
Problem=("1" "2" "3" "4" "5")
# ProblemSizeArray=("60")
PopulationSizeArray=("100" "300" "500" "1000" "2000")
# PopulationSizeArray=("30" "50")
SelectionPressureArray=("2" "5" "10" "15" "20")
# SelectionPressureArray=("15" "20")


# ell
for ((ell = 0; ell < 1; ell++))
do

  # Population Size
  for ((pop = 0; pop < 2; pop++))
  do

    # selection pressure
    for ((s = 0; s < 5; s++))
    do

  
      # # repeat times
      for ((times = 0; times < 10; times++))
      do
        
        
        # KMC
        # Problem
        for ((problem = 0; problem < 5; problem++))
        do
        
            ./RankingEDAsCEC -i ./LOP_instance/Cebe.lop.n100."${Problem[$problem]}"\
                              -o ./LOP_result/KMC/"${Problem[$problem]}"_"$((${PopulationSizeArray[$pop]}))"_"${SelectionPressureArray[$s]}"_KMC_"$times".txt\
                              -s "$times"\
                              -t LOP\
                              -m M -d C -v 0 -x 2\
                              -r "${SelectionPressureArray[$s]}"\
                              -z "$((${PopulationSizeArray[$pop]}))" &

        done
        # problem end


      done
      # repeat times end
      echo "end of pop_size = "$((${PopulationSizeArray[$pop]}))", Emax = "$((${ProblemSizeArray[$ell]} * ${ProblemSizeArray[$ell]} * 1000))", model = KMC, s = ${SelectionPressureArray[$s]}"
      wait

      
    done
    # selection pressure end



  done
   # Population Size end

done
# ell end