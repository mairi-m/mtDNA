import random
from mesa import Agent, Model
from mesa.datacollection import DataCollector
import matplotlib.pyplot as plt


class mtDNA(Agent):
    def __init__(self,unique_id,parent_id,r=0.01,d=0.01,status="wild-type"):
      super().__init__(unique_id, model)
        self.r=r
        self.d=d
        self.unique_id=unique_id
        self.parent_id=parent_id
        self.status=status
        
    def step(self):
      
      for mol_index in range(len(system_state)):
        molecule = system_state[mol_index]
      
        roll = random.random()
        if 0.0 < roll and roll < molecule.d:
            print("Degrading "+str(molecule.unique_id)+" at timepoint "+str(i)+"!")
          
            # Label molecule for removal
            self.model.molecules_to_remove.append(mol_index)

      
        
       
        
     

        
        
        
        
        
        
# Initialise cell
N0 = 20
system_state=[mtDNA(unique_id=x,parent_id=-1) for x in range(0,N0)]


class cell(Model):
  
  def __init__(self):
    self.schedule = RandomActivation(self)
    
    
    # Initialise cell
    N0 = 20
    system_state=[mtDNA(unique_id=x,parent_id=-1) for x in range(0,N0)]
    self.schedule.add(system_state)
    
    self.datacollector = DataCollector(
      model_reporters = {"Number of agents" : 
      agent_reporters = {"Status" : "status})
                                                
    
  def step(self):
    self.molecules_to_remove = []
    self.molecules_to_replicate = []                     
                         
    self.datacollector.collect(self)
    self.schedule.step()
                
    





# Simulation
current_id = N0 + 1
# Assume time step = 1
for i in range(0,10000): # like running the model 10000 times
    molecules_to_remove = []
    new_molecules = []
    
    for mol_ind in range(0,len(system_state)):
        molecule = system_state[mol_ind]
        
        roll = random.random()
        if 0.0 < roll and roll < molecule.d:
            print("Degrading "+str(molecule.unique_id)+" at timepoint "+str(i)+"!")
            # Label molecule for removal

            molecules_to_remove.append(mol_ind)
        elif molecule.d < roll and roll < molecule.r + molecule.d:
            print("Replicating "+str(molecule.unique_id)+" at timepoint "+str(i)+"!")
            # Label mother for removal
            molecules_to_remove.append(mol_ind)
            # Add two daughters
            for j in range(2):
                current_id += 1
                daughter = mtDNA(current_id, r=molecule.r, d=molecule.d, parent_id=molecule.unique_id)
                new_molecules.append(daughter)
                
    # Drop all molecules which have been labelled for removal
    system_state = [mol for i,mol in enumerate(system_state) if i not in molecules_to_remove]
    # Add in the newly formed daughter molecules
    system_state = system_state + new_molecules
    
