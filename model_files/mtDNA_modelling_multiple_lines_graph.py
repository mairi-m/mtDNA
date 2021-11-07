# Title: mtDNA modelling - Mutation load batch run (no MESA) - Copy number control
# Author: Mairi MacIain
# Date: 15th Oct 2021

'''
The base of the code was written by Conor and I am going to improve on it by adding
other things such as:

code:
Molecule generations
spatial effects
layering graphs


Analyse:
Record to file correlation with mutation loads and sim ending
mol_age effect
'''

# Imports
import random
import matplotlib.pyplot as plt

# Simulation function
def simulation(system_state, tmax = 5000):
    
    # Declare local variables
    system_states = []
    Nwt0 = len([mol for mol in system_state if mol.status=="wild-type"])
    Nmut0 = len([mol for mol in system_state if mol.status=="mutant"])
    current_id = Nwt0 + Nmut0
    times_to_record = list(range(0,tmax,10))

    # Assume time step = 1
    t = 0
    while t<tmax and len(system_state)>0:
        molecules_to_remove = []
        new_molecules = []
        if t in times_to_record:
            system_states.append(system_state)
            
        for mol_ind in range(0,len(system_state)):
            molecule = system_state[mol_ind]

            molecule.mol_age += 1
            
            roll = random.random()
            if roll <= molecule.d and len(system_state) > 30 and molecule.mol_age > 3:
                # Label molecule for removal
                molecules_to_remove.append(mol_ind)
                
            elif roll > molecule.d and roll <= (molecule.r + molecule.d) and len(system_state) < 90 and molecule.mol_age > 3:
                # Label mother for removal
                molecules_to_remove.append(mol_ind)
                
                # Add two daughters
                for j in range(2):
                    current_id += 1
                    daughter = mtDNA(current_id,parent_id=molecule.unique_id,mol_age=0,r=molecule.r,d=molecule.d,status=molecule.status)
                    new_molecules.append(daughter)
                    
        # Drop all molecules which have been labelled for removal
        system_state = [mol for i,mol in enumerate(system_state) if i not in molecules_to_remove]
        
        # Add in the newly formed daughter molecules
        system_state = system_state + new_molecules
        t = t+1
        
                    
    copy_numbers = [len(ss) for ss in system_states]
    # Remove the 0's from copy number array
    #copy_numbers = [no for no in copy_numbers if no != 0]

    mutants = [sum(mol.status=="mutant" for mol in ss) for ss in system_states]
    mutation_loads = [float(mutant)/float(copy_number) for mutant,copy_number in zip(mutants,copy_numbers)]

    sim_res = {"CN": copy_numbers, "ML": mutation_loads}
    return(sim_res)

# Function to make graph
def makePlot(copy_numbers, mutation_loads, tmax = 5000, filename = "mtDNA_pop.png", showplot=False):

    # Declare variables
    times_to_record_final = []
    times_to_record = range(0,tmax,10)

    for index in range(len(mutation_loads)):
        # To ensure the lists are of the same length
        len_t_r = len(times_to_record)
        len_m_l = len(mutation_loads[index])

        if len_t_r != len_m_l:

            print("Simulation not complete")

            difference = len_t_r - len_m_l    
            len_t_r = len_t_r - difference

            for e in range(len_t_r):
                times_to_record_final.append(times_to_record[e])

            times_to_record.append(times_to_record_final)
            times_to_record_final = []

        
    fig, (ax2, ax1) = plt.subplots(1, 2, figsize = (12, 6))
    fig.suptitle('Simulating mtDNA population dynamics')

    # Graph of mutation load/steps
    for pli in range(20):
        ax1.plot((times_to_record),mutation_loads[index])

    ax1.set_xlabel("Age (Steps)", fontsize = 15)
    ax1.set_ylabel("Mutation load", fontsize = 15)
    ax1.set_xticks(list(range(0,tmax,500)))
    ax1.set_xlim([0,tmax])

    # Graph of copy number/steps
    for plii in range(20):
        ax2.plot((times_to_record),copy_numbers[index])
        
    ax2.set_xlabel("Age (Steps)", fontsize = 15)
    ax2.set_ylabel("Total copy number", fontsize = 15)
    ax2.set_xticks(list(range(0,tmax,500)))
    ax2.set_xlim([0,tmax])

    if showplot:
        plt.show()
    else:
        plt.savefig(filename, dpi = 300, orientation = "landscape")

class mtDNA():
    def __init__(self,unique_id,parent_id,mol_age = 0,r=0.01,d=0.01,status="wild-type"):
        self.r=r
        self.d=d
        self.unique_id=unique_id
        self.parent_id=parent_id
        self.mol_age=mol_age
        self.status=status 

'''
MAIN CODE
'''
# Initialise cell
Nwt0 = 45
Nmut0 = 15
system_state = [mtDNA(unique_id=x,parent_id=-1,status="wild-type") for x in range(0,Nwt0)] + [mtDNA(unique_id=x,parent_id=-1,status="mutant") for x in range(0,Nmut0)]

# Declare plotting variables
times_to_record_final = []
times_to_record = range(0,tmax,10)
showplot = False


# Run simulation
results = [simulation(system_state) for i in range(0, 20)]
mutation_loads = results["ML"]
copy_numbers = results["CN"]

for index in range(len(mutation_loads)):
    # To ensure the lists are of the same length
     len_t_r = len(times_to_record)
     len_m_l = len(mutation_loads[index])

    if len_t_r != len_m_l:
        print("Simulation not complete")

        difference = len_t_r - len_m_l    
        len_t_r = len_t_r - difference

        for e in range(len_t_r):
            times_to_record_final.append(times_to_record[e])

        times_to_record.append(times_to_record_final)
        times_to_record_final = []

        
fig, (ax2, ax1) = plt.subplots(1, 2, figsize = (12, 6))
fig.suptitle('Simulating mtDNA population dynamics')

# Graph of mutation load/steps
for pli in range(20):
    ax1.plot((times_to_record),mutation_loads[index])

ax1.set_xlabel("Age (Steps)", fontsize = 15)
ax1.set_ylabel("Mutation load", fontsize = 15)
ax1.set_xticks(list(range(0,tmax,500)))
ax1.set_xlim([0,tmax])

# Graph of copy number/steps
for plii in range(20):
    ax2.plot((times_to_record),copy_numbers[index])
        
ax2.set_xlabel("Age (Steps)", fontsize = 15)
ax2.set_ylabel("Total copy number", fontsize = 15)
ax2.set_xticks(list(range(0,tmax,500)))
ax2.set_xlim([0,tmax])

if showplot:
    plt.show()
else:
    plt.savefig(filename, dpi = 300, orientation = "landscape")



