from mesa import Agent, Model
from mesa.time import RandomActivation
import numpy as np

class mtDNA(Agent):
    """An mtDNA molecule"""
    def __init__(self, unique_id, model, r, d, parent_id=0):
        super().__init__(unique_id, model)
        self.r = r
        self.d = d
    def step(self):
        roll = np.random.random()
        if 0.0 < roll and roll < self.d:
            self.model.remove_molecules.append(self)
        elif self.d < roll and roll < self.r + self.d:
            self.model.replicate_molecules.append(self)            
        

class cell(Model):
    def __init__(self, N, r=0.001, d=0.001):
        self.num_mols = N
        self.r = r
        self.d = d
        self.schedule = RandomActivation(self)
        for i in range(self.num_mols):
            molecule = mtDNA(i, self, r=self.r, d=self.d, parent_id=0)
            self.schedule.add(molecule)
            self.current_id = i

    def step(self):
        self.remove_molecules = []
        self.replicate_molecules = []
        self.schedule.step()
        for mol in self.replicate_molecules:
            print("Replicating "+str(mol.unique_id))
            for j in range(2):
                self.current_id += 1
                molecule = mtDNA(self.current_id, self, r=self.r, d=self.d, parent_id=mol.unique_id)
                self.schedule.add(molecule)
                self.num_mols += 1
        
        for mol in self.remove_molecules + self.replicate_molecules:
            print("Removing "+str(mol.unique_id))
            self.schedule.remove(mol)
            #self.remove_molecules.remove(mol)
            self.num_mols -= 1

mycell = cell(200)
# https://dmnfarrell.github.io/bioinformatics/abm-mesa-python
for i in range(1,100):
    mycell.step()

