import random
import matplotlib.pyplot as plt

class mtDNA():
    def __init__(self,unique_id,parent_id,r=0.01,d=0.01,status="wild-type"):
        self.r=r
        self.d=d
        self.unique_id=unique_id
        self.parent_id=parent_id
        self.status=status

# Initialise cell
Nwt0 = 15
Nmut0 = 5
tmax = 5000
system_state = [mtDNA(unique_id=x,parent_id=-1,status="wild-type") for x in range(0,Nwt0)] + [mtDNA(unique_id=x,parent_id=-1,status="mutant") for x in range(0,Nmut0)]

times_to_record = range(0,tmax,10)
system_states = []

# Simulation
current_id = Nwt0 + Nmut0
# Assume time step = 1
for t in range(0,tmax):
    if t in times_to_record:
        system_states.append(system_state)
    
    molecules_to_remove = []
    new_molecules = []
    
    for mol_ind in range(0,len(system_state)):
        molecule = system_state[mol_ind]
        
        roll = random.random()
        if 0.0 < roll and roll < molecule.d:
            #print("Degrading "+str(molecule.unique_id)+" at timepoint "+str(t)+"!")
            # Label molecule for removal

            molecules_to_remove.append(mol_ind)
        elif molecule.d < roll and roll < molecule.r + molecule.d:
            #print("Replicating "+str(molecule.unique_id)+" at timepoint "+str(t)+"!")
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
mutants = [sum(mol.status=="mutant" for mol in ss) for ss in system_states]
mutation_loads = [float(mutant)/float(copy_number) for mutant,copy_number in zip(mutants,copy_numbers)]

plt.plot(times_to_record,mutation_loads)
plt.show()

plt.plot(times_to_record,copy_numbers)
plt.show()
            
