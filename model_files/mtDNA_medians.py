# Title: Modelling mtDNA population dynamics
# Author: Mairi MacIain
# Date: 15th Oct 2021


# Imports
import random
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import csv



''' Simulation function'''


def simulation(system_state, times_to_record, tmax = 5000, CN_upper = 90, CN_lower = 30):
    
    # Declare local variables
    system_states = []
    Nwt0 = len([mol for mol in system_state if mol.status=="wild-type"])
    Nmut0 = len([mol for mol in system_state if mol.status=="mutant"])
    current_id = Nwt0 + Nmut0
    #times_to_record = list(range(0,tmax,10))

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
            rand_age = random.randint(2,5)
            
            roll = random.random()
            if roll <= molecule.d and len(system_state) > CN_lower and molecule.mol_age > rand_age:
                # Label molecule for removal
                molecules_to_remove.append(mol_ind)
                
            elif roll > molecule.d and roll <= (molecule.r + molecule.d) and len(system_state) < CN_upper and molecule.mol_age > rand_age:
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


    mutants = [sum(mol.status=="mutant" for mol in ss) for ss in system_states]
    mutation_loads = [float(mutant)/float(copy_number) for mutant,copy_number in zip(mutants,copy_numbers)]

    sim_res = {"CN": copy_numbers, "ML": mutation_loads}
    return(sim_res)



''' Function to calculate percentiles'''

def percentiles(results, var = "ML"):
    values = []
    d50 = []


    for i in range(len(results[0]["ML"])): # All have the same length because of CN control
        for j in range(len(results)):
            values.append(results[j][var][i])

        d50.append(np.percentile(values, 50))



    return(d50)
        
    
## Simulate median mutation load dynamics for lots of initial conditions
def sim_med(Nwt0, Nmut0,  times_to_record, runs=35, tmax = 5000):

    ML_distributions = []

    for i in range(len(Nwt0)):
        system_state = [mtDNA(unique_id=x,parent_id=-1,status="wild-type") for x in range(0,Nwt0[i])] + [mtDNA(unique_id=x,parent_id=-1,status="mutant") for x in range(0,Nmut0[i])]
        results = [simulation(system_state, times_to_record, tmax, CN_upper = ((Nwt0[i]+Nmut0[i])*1.5), CN_lower =((Nwt0[i]+Nmut0[i])*0.5)) for j in range(0, runs)]
        ML_distributions.append(percentiles(results, "ML"))

    return(ML_distributions)


def plot_med(times_to_record, sim_meds, Nwt0, filename):

    times_to_record_y = [tr/365.0 for tr in times_to_record]

    for p,sm in enumerate(sim_meds):
        plt.plot(times_to_record_y, sm, color = "red")
        IC = Nmut0[p] / (Nmut0[p] + Nwt0[p])
        plt.axhline(y = IC)
    
    plt.title("Medians of mutation loads")
    plt.xlabel("Age (years)", fontsize = 15)
    plt.ylabel("Mutation load", fontsize = 15)
    plt.xticks(list(range(0,round(tmax/365.0),5)))
    plt.xlim([0,tmax/365])
    plt.ylim([0.0, 1.0])

    plt.savefig(filename, dpi = 300, orientation = "landscape")



''' END OF FUNCTIONS'''




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
tmax = 30000
runs = 300
times_to_record = range(0,tmax,10)
Nwt0 = [200, 180, 160, 140, 120, 100, 80, 60, 40, 20, 0]
Nmut0 = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200]

dirname = "Median_Analysis"
Path(dirname).mkdir(parents = True, exist_ok = True)


sim_meds = sim_med(Nwt0,Nmut0,times_to_record,runs,tmax)

# Save file to hard disk
with open("sim_meds.csv", "w",newline="") as f:
    wr = csv.writer(f)
    wr.writerows(sim_meds)

# Read them back in again (everything from here down can be moved to a separate file)
with open("sim_meds.csv", "r") as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(read_obj, quoting=csv.QUOTE_NONNUMERIC)
    # Pass reader object to list() to get a list of lists
    sim_meds = list(csv_reader)


plot_med(times_to_record, sim_meds, Nwt0, filename = os.path.join(dirname,"medians(2).png"))















    
