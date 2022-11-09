


# problem: t65b11xx, t65d11xx, t65f11xx, satbu1, satbu2, satbu3, N-t65b11xx_150, N-t65d11xx_150, N-t65f11xx_150, N-t65l11xx_150
# model: EO, ET, NO, NT


population_size=("40" "80" "160" "320" "640")

instance=("t65b11xx" "t65d11xx" "t65f11xx" "satbu1" "satbu2" "satbu3" "N-t65b11xx_150" "N-t65d11xx_150" "N-t65f11xx_150" "N-t65l11xx_150")
ell=("44" "44" "44" "60" "60" "60" "150" "150" "150" "150")

# problem num
for ((prob = 0; prob < 10; prob++))
do

    # pop size
    for ((pop = 0; pop < 5; pop++))
    do


        # LOP
        # repeat times
        for ((times = 0; times < 10; times++))
        do

            echo "problem = "${instance[$prob]}", pop_size = "${population_size[$pop]}", times = ${times}, Emax = "$((${ell[$prob]} * ${ell[$prob]} * 1000))" "

            ./RankingEDAsCEC -i ./LOP_instance/"${instance[$prob]}"\
                              -o ./LOP_result/KMC/"${instance[$prob]}"_"${population_size[$pop]}"_KMC_"$times".txt\
                              -s "$times"\
                              -t LOP\
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
