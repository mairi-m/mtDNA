# Title: Modelling mtDNA population dynamics
# Author: Mairi MacIain
# Date: 15th Oct 2021


# Imports
import random
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os



''' Simulation function'''


def simulation(system_state, tmax = 5000, CN_upper = 90, CN_lower = 30):
    
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
    d5 = []
    d25 = []
    d50 = []
    d75 = []
    d95 = []

    for i in range(len(results[0]["ML"])): # All have the same length because of CN control
        for j in range(len(results)):
            values.append(results[j][var][i])

        d5.append(np.percentile(values, 5))
        d25.append(np.percentile(values, 25))
        d50.append(np.percentile(values, 50))
        d75.append(np.percentile(values, 75))
        d95.append(np.percentile(values, 95))

    all_distributions = {"d5": d5, "d25": d25, "d50": d50, "d75": d75, "d95": d95}
    return(all_distributions)
        
    


''' Function to make plots'''


def makePlot(results, ML_distributions, CN_distributions, Ntot0, runs, tmax = 5000, filename = "Multiple lines plots/mtDNA_pop.png", showplot = False, showboth = False):

    # Declare variables
    times_to_record_final = []
    times_to_record = range(0,tmax,10)



    times_to_record_y = [tr/365.0 for tr in times_to_record]

    
        
    fig, (ax2, ax1) = plt.subplots(1, 2, figsize = (12, 6))
    fig.suptitle('Simulating mtDNA population dynamics')

    # Graph of mutation load/steps
    for plot in range(runs):
        ax1.plot((times_to_record_y),results[plot]["ML"], color = "black", alpha = 0.08)
    ax1.plot((times_to_record_y), ML_distributions["d5"], color = "r", ls = ":")
    ax1.plot((times_to_record_y), ML_distributions["d25"], color = "r", ls = "--")
    ax1.plot((times_to_record_y), ML_distributions["d50"], color = "r", linewidth = "3")
    ax1.plot((times_to_record_y), ML_distributions["d75"], color = "r", ls = "--")
    ax1.plot((times_to_record_y), ML_distributions["d95"], color = "r", ls = ":")
    
    ax1.set_xlabel("Age (years)", fontsize = 15)
    ax1.set_ylabel("Mutation load", fontsize = 15)
    ax1.set_xticks(list(range(0,round(tmax/365.0),1)))
    ax1.set_xlim([0,tmax/365])
    ax1.set_ylim([0.0, 1.0])

    # Graph of copy number/steps
    for plot in range(runs):
        ax2.plot((times_to_record_y),results[plot]["CN"], color = "black", alpha = 0.08)
    ax2.plot((times_to_record_y), CN_distributions["d5"], color = "r", ls = ":")
    ax2.plot((times_to_record_y), CN_distributions["d25"], color = "r", ls = "--")
    ax2.plot((times_to_record_y), CN_distributions["d50"], color = "r", linewidth = "3")
    ax2.plot((times_to_record_y), CN_distributions["d75"], color = "r", ls = "--")
    ax2.plot((times_to_record_y), CN_distributions["d95"], color = "r", ls = ":")
    
    ax2.set_xlabel("Age (years)", fontsize = 15)
    ax2.set_ylabel("Total copy number", fontsize = 15)
    ax2.set_xticks(list(range(0,round(tmax/365.0),1)))
    ax2.set_xlim([0,tmax/365])
    ax2.set_ylim([0.0,((Ntot0 *1.5) + 50.0)])

    if showplot:
        plt.show()
    elif showboth:
        plt.savefig(filename, dpi = 300, orientation = "landscape")  
        plt.show()

    else:
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
tmax = 5000
runs = 35
times_to_record = range(0,tmax,10)
Nwt0 = 50
Nmut0 = 150
Ntot0 = Nwt0 + Nmut0
system_state = [mtDNA(unique_id=x,parent_id=-1,status="wild-type") for x in range(0,Nwt0)] + [mtDNA(unique_id=x,parent_id=-1,status="mutant") for x in range(0,Nmut0)]


# Run simulation
results = [simulation(system_state, CN_upper = ((Nwt0+Nmut0)*1.5), CN_lower =((Nwt0+Nmut0)*0.5)) for i in range(0, runs)]

dirname = "Analysis"
Path(dirname).mkdir(parents = True, exist_ok = True)

# Get distributions
ML_distributions = percentiles(results, "ML")
CN_distributions = percentiles(results, "CN")

# Make plots
makePlot(results, ML_distributions, CN_distributions, Ntot0, runs, filename = os.path.join(dirname,"mtDNA_pop.png"))


