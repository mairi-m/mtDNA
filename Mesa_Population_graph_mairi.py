import random
from mesa import Agent, Model
from mesa.time import RandomActivation
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
            print("Degrading "+str(molecule.unique_id))
            # Label molecule for removal
            self.model.molecules_to_remove.append(mol_index)
            
         elif molecule.d < roll and roll < molecule.r + molecule.d:
            print("Replicating "+str(molecule.unique_id))
            # Label mother for removal
            self.model.molecules_to_remove.append(mol_index)
            self.model.mother_mol_uniqueID(molecule.unique_id) # save mother id to give to daughter in model class

      
        
       


class cell(Model):
  
  def __init__(self):
    self.schedule = RandomActivation(self)
    
    
    # Initialise cell
    N0 = 20
    system_state=[mtDNA(unique_id=x,parent_id=-1) for x in range(0,N0)]
    self.schedule.add(system_state)
    
    self.datacollector = DataCollector(
      model_reporters = {"Number of agents" : len(system_state)}
      agent_reporters = {"Status" : "status})
                                                
    
  def step(self):    
    current_id = N0
    self.molecules_to_remove = []
    self.mother_mol_uniqueID= []                     
    self.new_molecules = []

    for index in range(len(self.mother_mol_uniqueID):
                       for j in range(2):
                       current_id += 1
                       daughter = mtDNA(current_id, r = 0.01, d = 0.01, parent_id = self.mother_mol_uniqueID[index]
                       self.new_molecules.append(daughter)
                       
     
                                        
                
    # Drop all molecules which have been labelled for removal
    system_state = [mol for i,mol in enumerate(system_state) if i not in molecules_to_remove]
    # Add in the newly formed daughter molecules
    system_state = system_state + new_molecules
                         
                         
                     
    self.datacollector.collect(self)
    self.schedule.step()
                                        
                                        
   model = cell():
  for index in range(100):
           model.step()
                                        
  no_of_mtDNA = model.datacollector.get_model_vars_dataframe()
no_of_mtDNA.plot()
plt.show()
                
    


    
