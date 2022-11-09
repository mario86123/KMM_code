
# def read_avg_fitness(problem, model, pop_size):
    
#     summation = 0
    
#     for t in range(10):
#         f = open(f'./{model}/{problem}_{pop_size}_{model}_{t}.txt', 'r')
#         s = f.read().splitlines()

#         best_fitness = float(s[0].split(": ")[1])
#         # best_fitness_NFE = float(s[2].split(": ")[1])

#         summation += best_fitness


#     return summation / 10; 

def read_arpd(problem, selection_pressure, pop_size, opt):
    
    summation = 0
    
    for t in range(10):
        # print(f'./{model}/{problem}_{pop_size * ell}_{model}_{t}.txt')
        f = open(f'./KMC/{problem}_{pop_size}_{selection_pressure}_KMC_{t}.txt', 'r')
        s = f.read().splitlines()

        best_fitness = float(s[0].split(": ")[1])
        # best_fitness_NFE = float(s[2].split(": ")[1])
        
        rpd = ((float(-best_fitness)) - (-opt)) / (-opt) * 100

        # print(rpd)
        summation += rpd


    return summation / 10; 



# def read_avg_NFE(problem, model, pop_size):

#     summation = 0
    
#     for t in range(10):
#         f = open(f'./{model}/{problem}_{pop_size}_{model}_{t}.txt', 'r')
#         s = f.read().splitlines()

#         # best_fitness = float(s[0].split(": ")[1])
#         best_fitness_NFE = float(s[2].split(": ")[1])

#         summation += best_fitness_NFE


#     return summation / 10; 



problem_instance_lst = ["100_5", "100_10", "100_20"]
# ell_lst = [20, 50]
opt_lst = [-253713, -299431, -367267] # according to EHBSA paper




pop_size_lst = [100, 300, 500, 1000, 2000]
# pop_size_lst = [1, 3, 5, 10, 20, 40, 80, 160]
# other_pop_size_lst = [5, 10, 20]

# model_lst = ["MST", "GMC", "ET5", "NO",]
# model_lst = ["MST", "GMC", "MU", "ET5", "NO",]
selection_pressure_lst = [2, 5, 10, 15, 20]


print()
print(" [[[ ARPD lower, better ]]] ")
print()
print("              |  s=2  |  s=5  |  s=10  |  s=15  |  s=20  |")
# print("            |  MST  |  GMC  |  ET5  |   NO  |")
# print("              |  MR  | MST | MSTME|  GMC |  MU  |  ET5 |")

for ell_count in range(len(problem_instance_lst)):
    
    print("-----------------------------------------")

    for pop_size_count in range(len(pop_size_lst)):

        print(f'{problem_instance_lst[ell_count]} |', end=" ")
        print(f'{(pop_size_lst[pop_size_count]):5} |', end=" ")

        for selection_pressure_count in range(len(selection_pressure_lst)):
                print(f'{read_arpd(problem_instance_lst[ell_count] , selection_pressure_lst[selection_pressure_count], pop_size_lst[pop_size_count], opt_lst[ell_count]):5.1f} |', end=" ")
            
            # else:
                # print("      |", end=" ")









        print()

print("------------------------------------------------------------------------------------")
# read_avg_fitness(problem, model, ell, pop_size):


