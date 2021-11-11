# Title: mtDNA modelling - Mutation load batch run (no MESA) - Copy number control
# Author: Mairi MacIain
# Date: 15th Oct 2021

'''
The base of the code was written by Conor and I am going to improve on it by adding
other things such as:



Analyse:
Record to file correlation with mutation loads and sim ending
mol_age effect
'''

# Imports
import random
import numpy
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
    # Remove the 0's from copy number array
    #copy_numbers = [no for no in copy_numbers if no != 0]

    mutants = [sum(mol.status=="mutant" for mol in ss) for ss in system_states]
    mutation_loads = [float(mutant)/float(copy_number) for mutant,copy_number in zip(mutants,copy_numbers)]

    sim_res = {"CN": copy_numbers, "ML": mutation_loads}
    return(sim_res)



''' Function to make plots'''


def makePlot(results, ML_distributions, CN_distributions, tmax = 5000, filename = "Multiple lines plots/mtDNA_pop_20sims.png", showplot = False):

    # Declare variables
    times_to_record_final = []
    times_to_record = range(0,tmax,10)


        
    fig, (ax2, ax1) = plt.subplots(1, 2, figsize = (12, 6))
    fig.suptitle('Simulating mtDNA population dynamics')

    # Graph of mutation load/steps
    for plot in range(20):
        ax1.plot((times_to_record),results[plot]["ML"], color = "black", alpha = 0.3)
    ax1.plot((times_to_record), ML_distributions["M5"], color = "g", ls = "--")
    ax1.plot((times_to_record), ML_distributions["M25"], color = "g", ls = "--")
    ax1.plot((times_to_record), ML_distributions["M50"], color = "g", linewidth = "3")
    ax1.plot((times_to_record), ML_distributions["M75"], color = "g", ls = "--")
    ax1.plot((times_to_record), ML_distributions["M95"], color = "g", ls = "--")
    
    ax1.set_xlabel("Age (Steps)", fontsize = 15)
    ax1.set_ylabel("Mutation load", fontsize = 15)
    ax1.set_xticks(list(range(0,tmax,500)))
    ax1.set_xlim([0,tmax])
    ax1.set_ylim([0.0, 1.0])

    # Graph of copy number/steps
    for plot in range(20):
        ax2.plot((times_to_record),results[plot]["CN"], color = "black", alpha = 0.3)
    ax2.plot((times_to_record), CN_distributions["C5"], color = "g", ls = "--")
    ax2.plot((times_to_record), CN_distributions["C25"], color = "g", ls = "--")
    ax2.plot((times_to_record), CN_distributions["C50"], color = "g", linewidth = "3")
    ax2.plot((times_to_record), CN_distributions["C75"], color = "g", ls = "--")
    ax2.plot((times_to_record), CN_distributions["C95"], color = "g", ls = "--")
    
    ax2.set_xlabel("Age (Steps)", fontsize = 15)
    ax2.set_ylabel("Total copy number", fontsize = 15)
    ax2.set_xticks(list(range(0,tmax,500)))
    ax2.set_xlim([0,tmax])
    ax2.set_ylim([0.0,100.0])

    if showplot:
        plt.show()
    else:
        plt.savefig(filename, dpi = 300, orientation = "landscape")


def ML_percentiles(results, times_to_record):

    mutation_loads = []
    Mdistribution5 = []
    Mdistribution25 = []
    Mdistribution50 = []
    Mdistribution75 = []
    Mdistribution95 = []

    for i in range(len(times_to_record)):

        for j in range(20):
            mutation_loads.append(results[j]["ML"][i])

        dist5 = numpy.percentile(mutation_loads, 5)
        dist25 = numpy.percentile(mutation_loads, 25)
        dist50 = numpy.percentile(mutation_loads, 50)
        dist75 = numpy.percentile(mutation_loads, 75)
        dist95 = numpy.percentile(mutation_loads, 95)

        Mdistribution5.append(dist5)
        Mdistribution25.append(dist25)
        Mdistribution50.append(dist50)
        Mdistribution75.append(dist75)
        Mdistribution95.append(dist95)

    all_distributions = {"M5": Mdistribution5, "M25": Mdistribution25, "M50": Mdistribution50, "M75": Mdistribution75, "M95": Mdistribution95,}
    return(all_distributions)


def CN_percentiles(results, times_to_record):

    copy_number = []
    Cdistribution5 = []
    Cdistribution25 = []
    Cdistribution50 = []
    Cdistribution75 = []
    Cdistribution95 = []

    for i in range(len(times_to_record)):

        for j in range(20):
            copy_number.append(results[j]["CN"][i])

        dist5 = numpy.percentile(copy_number, 5)
        dist25 = numpy.percentile(copy_number, 25)
        dist50 = numpy.percentile(copy_number, 50)
        dist75 = numpy.percentile(copy_number, 75)
        dist95 = numpy.percentile(copy_number, 95)

        Cdistribution5.append(dist5)
        Cdistribution25.append(dist25)
        Cdistribution50.append(dist50)
        Cdistribution75.append(dist75)
        Cdistribution95.append(dist95)

    all_distributions = {"C5": Cdistribution5, "C25": Cdistribution25, "C50": Cdistribution50, "C75": Cdistribution75, "C95": Cdistribution95,}
    return(all_distributions)


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
times_to_record = range(0,tmax,10)
Nwt0 = 45
Nmut0 = 15
system_state = [mtDNA(unique_id=x,parent_id=-1,status="wild-type") for x in range(0,Nwt0)] + [mtDNA(unique_id=x,parent_id=-1,status="mutant") for x in range(0,Nmut0)]


# Run simulation
results = [simulation(system_state) for i in range(0, 20)]

dirname = "Multiple lines plots"
Path(dirname).mkdir(parents = True, exist_ok = True)

# Get distributions
ML_distributions = ML_percentiles(results, times_to_record)
CN_distributions = CN_percentiles(results, times_to_record)

# Make plots
makePlot(results, ML_distributions, CN_distributions, filename = os.path.join(dirname,"mtDNA_pop_20sims.png"))


