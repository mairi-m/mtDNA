# Title: mtDNA population graph - no mutations
# Author: Mairi MacIain
# Date: 29th Sept 2021

''' Imports '''
import random
import matplotlib.pyplot as plt
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.datacollection import DataCollector


''' Create agent ''' # These multi-line string comments are only useful in the context of first line of function
class mtDNA(Agent):
    def __init__(self, unique_id, model, parent_id, r = 0.01, d = 0.01, status = "wild-type"):
        super().__init__(unique_id, model)
        self.r = r
        self.d = d
        self.unique_id = unique_id
        self.parent_id = parent_id
        self.status = status

    def step(self): 

        ''' Iterate through the agents '''
        for mol_index in range (len(system_state)):
            molecule = system_state[mol_index]
            
            ''' If the molecule is to be degraded add it to the list for removal''' # These multi-line string comments are only useful in the context of first line of function
            roll = random.random()
            if 0.0 < roll and roll < molecule.d:
                print("Degrading " + str(molecule.unique_id))
                self.model.molecules_to_remove.append(mol_index)
                
                ''' If the molecule is to be replicated add it to the list for removal and add its unique id
to a list so the amount of molecules to be replicated, and their parent Ids are known '''
            elif molecule.d < roll and roll < molecule.r + molecule.d:
                print("Replicating" + str(molecule.unique_id))
                self.model.molecules_to_remove.append(mol_index)
                self.model.mother_mol_uniqueID.append(molecule.unique_id)



''' Create the model'''
class cell(Model):

    def __init__(self, N0 = 20):
        self.schedule = RandomActivation(self)
        
        ''' Initialise the cell '''
        system_state = [mtDNA(x, self, parent_id = -1) for x in range(0, N0)] # This is not compatible with MESA (this is from the simulation engine we wrote from scratch)
        # Need to be able to add "self" here as an argument to mtDNA, in order to make the connection between agent and model.  
        self.schedule.add(system_state) # This is where the incompatibility bites.  self.schedule.add expects a single mtDNA agent, whereas you are trying to pass a list
        
        ''' Create the data collector'''
        self.datacollector = DataCollector(
            model_reporters = {"Number of agents" : len(system_state)},
            agent_reporters = {"Status" : "status"})

    def step(self):
        ''' Declare variables'''
        current_id = N0
        self.molecules_to_remove = []
        self.mother_mol_uniqueID = []
        self.new_molecules = []
        
        ''' Create new daughter molecules with new ids and their parents id'''
        for index in range(len(self.mother_mol_uniqueID)):
            for j in range(2):
                current_id += 1
                daughter = mtDNA(self,current_id, parent_id = self.mother_mol_uniqueID[index], r = 0.01, d = 0.01, status = "wild-type")
                self.new_molecules.append(daughter)
                
        ''' Remove the degraded and mother molecules'''
        system_state = [mol for i, mol in enumerate(system_state) if i not in molecules_to_remove] # You have added the daughters but haven't removed the mother molecules
        system_state = system_state + new_molecules


        self.datacollector.collect(self)
        self.schedule.step()


model = cell()
for index in range(100):
        model.step()


mtDNA_pop_graph = model.datacollector.get_model_vars_dataframe()
mtDNA_pop_graph.plot()
plt.show()

