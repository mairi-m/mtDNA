import random

class mtDNA():
    def __init__(self,unique_id,parent_id,r=0.01,d=0.01,status="wild-type"):
        self.r=r
        self.d=d
        self.unique_id=unique_id
        self.parent_id=parent_id
        self.status=status

# Initialise cell
system_state=[mtDNA(unique_id=x,parent_id=-1) for x in range(0,20)]

# Simulation
# Assume time step = 1
for i in range(0,100):
    for molecule in system_state:
        molecules_to_remove
        roll = random.random()
        if 0.0 < roll and roll < cell.d:
            # Remove molecule from system_state
            system_state.remove(molecule)
        elif molecule.d < roll and roll < molecule.r + molecule.d:
            # Replicate molecule
            for j in range(2):
                current_id += 1
                daughter = mtDNA(current_id, r=molecule.r, d=molecule.d, parent_id=molecule.unique_id)
                system_state.append(daughter)
                
            
            
            
