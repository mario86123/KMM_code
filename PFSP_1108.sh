


# problem: 50_5_1, 50_5_2, 50_10_1, 50_10_2, 50_20_1, 50_20_2, 100_5_1, 100_5_2, 100_10_1, 100_10_2, 100_20_1, 100_20_2

# model: EO, ET, NO, NT


# pop size:
# WO: 40, 80, 160, 320, 640
# WT: 10, 20, 40, 80, 160


population_size=("40" "80" "160" "320" "640")

instance=("50_5_1" "50_5_2" "50_10_1" "50_10_2" "50_20_1" "50_20_2" "100_5_1" "100_5_2" "100_10_1" "100_10_2" "100_20_1" "100_20_2")
ell=("50" "50" "50" "50" "50" "50" "100" "100" "100" "100" "100" "100")


# problem num
for ((prob = 0; prob < 12; prob++))
do

    # pop size
    for ((pop = 0; pop < 5; pop++))
    do


        # PFSP
        # repeat times
        for ((times = 0; times < 10; times++))
        do

            echo "problem = "${instance[$prob]}", pop_size = "${population_size[$pop]}", times = ${times}, Emax = "$((${ell[$prob]} * ${ell[$prob]} * 1000))" "


            ./RankingEDAsCEC -i ./PFSP_instance/tai"${instance[$prob]}".txt\
                              -o ./PFSP_result/KMC/"${instance[$prob]}"_"${population_size[$pop]}"_KMC_"$times".txt\
                              -s "$times"\
                              -t PFSP\
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
