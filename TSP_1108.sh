


# problem: gr24, gr48, eil51, berlin52, eil76, pr76, gr96, rat99, eil101, pr107
# model: EO, ET, NO, NT


# pop size:
# WO: 40, 80, 160, 320, 640
# WT: 10, 20, 40, 80, 160


population_size=("40" "80" "160" "320" "640")

instance=("gr24" "gr48" "eil51" "berlin52" "eil76" "pr76" "gr96" "rat99" "eil101" "pr107")
ell=("24" "48" "51" "52" "76" "76" "96" "99" "101" "107")


# problem num
for ((prob = 0; prob < 10; prob++))
do

    # pop size
    for ((pop = 0; pop < 5; pop++))
    do


        # TSP
        # repeat times
        for ((times = 0; times < 10; times++))
        do

            echo "problem = "${instance[$prob]}", pop_size = "${population_size[$pop]}", times = ${times}, Emax = "$((${ell[$prob]} * ${ell[$prob]} * 1000))" "

            ./RankingEDAsCEC -i ./TSP_instance/"${instance[$prob]}".tsp\
                              -o ./TSP_result/KMC/"${instance[$prob]}"_"${population_size[$pop]}"_KMC_"$times".txt\
                              -s "$times"\
                              -t TSP\
                              -m M -d C -v 0 -x 2\
                              -r 10\
                              -z "${population_size[$pop]}" &


        done
        # repeat times end
        wait

    done
    wait
    # pop size end


done
wait
# problem num end
