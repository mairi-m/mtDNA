# Title: mtDNA modelling - Mutation load batch run (no MESA)
# Author: Mairi MacIain
# Date: 15th Oct 2021

'''
The base of the code was written by Conor and I am going to improve on it by adding
other things such as:
Batch runs
Changing stats for mutants
Allowing wild type molecules to mutate
Molecule ages
Molecule generations


make tmax a changable variable
change fonts
'''

# Imports
import random
import matplotlib.pyplot as plt



# Simulation function
def simulation(system_state, Nwt0, Nmut0, tmax = 5000):
    
    # Declare local variables
    system_states = []
    molecules_to_remove = []
    new_molecules = []
    current_id = Nwt0 + Nmut0
    times_to_record = range(0,tmax,10)

    
    # Assume time step = 1
    for t in range(0,tmax):
        if t in times_to_record:
            system_states.append(system_state)
        

        
        for mol_ind in range(0,len(system_state)):
            molecule = system_state[mol_ind]
            
            roll = random.random()
            if 0.0 < roll and roll < molecule.d:
                # Label molecule for removal
                molecules_to_remove.append(mol_ind)
                
            elif molecule.d < roll and roll < molecule.r + molecule.d:
                # Label mother for removal
                molecules_to_remove.append(mol_ind)
                
                # Add two daughters
                for j in range(2):
                    current_id += 1
                    daughter = mtDNA(current_id, r=molecule.r, d=molecule.d, parent_id=molecule.unique_id, status=molecule.status)
                    new_molecules.append(daughter)
                    
        # Drop all molecules which have been labelled for removal
        system_state = [mol for i,mol in enumerate(system_state) if i not in molecules_to_remove]
        
        # Add in the newly formed daughter molecules
        system_state = system_state + new_molecules
        
                    
    copy_numbers = [len(ss) for ss in system_states]
    # Remove the 0's from copy number array
    copy_numbers = [no for no in copy_numbers if no != 0]

    mutants = [sum(mol.status=="mutant" for mol in ss) for ss in system_states]
    mutation_loads = [float(mutant)/float(copy_number) for mutant,copy_number in zip(mutants,copy_numbers)]

    
    sim_res = {"CN": copy_numbers, "ML": mutation_loads}
    return(sim_res)


# Function to make graph
def makePlot(copy_numbers, mutation_loads, tmax = 5000, filename = "mtDNA_pop.png"):

    # Declare variables
    times_to_record_final = []
    times_to_record = range(0,tmax,10)

    # To ensure the lists are of the same length
    len_t_r = len(times_to_record)
    len_m_l = len(mutation_loads)

    if len_t_r != len_m_l:

        print("Simulation not complete")

        difference = len_t_r - len_m_l    
        len_t_r = len_t_r - difference

        for e in range(len_t_r):
            times_to_record_final.append(times_to_record[e])  # Times to record needs emptied

        times_to_record = times_to_record_final
        times_to_record_final = []

        
    fig, (ax2, ax1) = plt.subplots(1, 2, figsize = (12, 6))
    fig.suptitle('Simulating mtDNA population dynamics')

    # Graph of mutation load/steps
    ax1.plot(times_to_record,mutation_loads)
    ax1.set_xlabel("Age (steps)", fontsize = 15)
    ax1.set_ylabel("Mutation load", fontsize = 15)
    ax1.set_xticks(list(range(0,5500,500)))
    ax1.set_xlim([0,5500])

    # Graph of copy number/steps
    ax2.plot(times_to_record,copy_numbers)
    ax2.set_xlabel("Age (steps)", fontsize = 15)
    ax2.set_ylabel("Total copy number", fontsize = 15)
    ax2.set_xticks(list(range(0,5500,500)))
    ax2.set_xlim([0,5500])



    plt.show()
    #plt.savefig(filename, dpi = 300, orientation = "landscape")

    



'''
MAIN CODE
'''
# Initialise cell
class mtDNA():
    def __init__(self,unique_id,parent_id,r=0.01,d=0.01,status="wild-type"):
        self.r=r
        self.d=d
        self.unique_id=unique_id
        self.parent_id=parent_id
        self.status=status

Nwt0 = 15
Nmut0 = 5
system_state = [mtDNA(unique_id=x,parent_id=-1,status="wild-type") for x in range(0,Nwt0)] + [mtDNA(unique_id=x,parent_id=-1,status="mutant") for x in range(0,Nmut0)]
    


# Run simulation 

sim_res = simulation(system_state, Nwt0, Nmut0)
mutation_loads = sim_res["ML"]
copy_numbers = sim_res["CN"]

results = [simulation(system_state) for i in range(0, 3)]

for i,sim_res in enumerate(results):
    mutation_loads = sim_res["ML"]
    copy_numbers = sim_res["CN"]
    makePlot(copy_numbers, mutation_loads, filename = "mtDNA_" + str(i) + ".png")
















