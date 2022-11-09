


# problem: bur26a, bur26b, tai40a, tai40b, tai60a, tai60b, tai80a, tai80b, tai100a, tai100b
# model: EO, ET, NO, NT


# pop size:
# WO: 40, 80, 160, 320, 640
# WT: 10, 20, 40, 80, 160


population_size=("40" "80" "160" "320" "640")

instance=("bur26a" "bur26b" "tai40a" "tai40b" "tai60a" "tai60b" "tai80a" "tai80b" "tai100a" "tai100b")
ell=("26" "26" "40" "40" "60" "60" "80" "80" "100" "100")


# problem num
for ((prob = 0; prob < 10; prob++))
do

    # pop size
    for ((pop = 0; pop < 5; pop++))
    do


        # QAP
        # repeat times
        for ((times = 0; times < 10; times++))
        do

            echo "problem = "${instance[$prob]}", pop_size = "${population_size[$pop]}", times = ${times}, Emax = "$((${ell[$prob]} * ${ell[$prob]} * 1000))" "

            ./RankingEDAsCEC -i ./QAP_instance/"${instance[$prob]}".dat.dat\
                              -o ./QAP_result/KMC/"${instance[$prob]}"_"${population_size[$pop]}"_KMC_"$times".txt\
                              -s "$times"\
                              -t QAP\
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
