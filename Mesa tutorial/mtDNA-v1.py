# Title: mtDNA v1



# Imports
import time
import enum
import numpy as np
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector



class Cell(Model):

        def __init__(self, N = 30, width = 20, height = 20, probMutant = 0.01):

                self.num_agents = N
                self.cellSize = MultiGrid(width, height, True)
                self.probabilityMutant = probMutant
                self.schedule = RandomActivation(self)

                # Create agents

                for index in range(self.num_agents):
                        a = mtDNA(self, index)
                        self.schedule.add(a)

                        # Add agent to a place in the grid (the cell)

                        x = self.random.randrange(self.cellSize.width)
                        y = self.random.randrange(self.cellSize.height)
                        self.cellSize.place_agent(a, (x,y))

                        # Make some mutant
                        mutant = np.random.randrange(0, 1)
                        if mutant < self.probabilityMutant:
                                a.status = Status.Mutant

                self.datacollector = DataCollector(
                        agent_reporters = {"Status": "status"}
                )

        def step(self):
                self.datacollector.collect(self)
                self.schedule.step()

class Status(enum.IntEnum):
        Mutant = 0
        Wild = 1
        Dead = 2

class mtDNA(Agent):

        def __init__(self, unique_id, model):
                super().__init__(unique_id, model)
                self.age = 0
                self.status = Status.Wild

        def move(self):
                possible_steps = self.model.grid.get_neighborhood(
                        self.pos,
                        moore = True,
                        include_center = True
                )

                new_position = self.random.choice(possible_steps)
                self.model.grid.move_agent(self, new_position)


        # Check age and replicate / mutate / remove

        def age_increase(self):

                self.age = self.age + 1

# sort out amount of agents because if  the keep creating every step then replicating it is not accurate

        def replicate(self):

                new_molecule_id = 0
                # Find an age for the molecule to replicate
                rep_age = self.random.randrange(20,25)

                # Check if they are of replication age
                if self.age > rep_age:
                        self.status = Status.Dead
                        for i in range(2):
                                a = mtDNA(self, new_molecule_id)
                                # Make some mutant
                                mutant = np.random.randrange(0, 1)
                                if mutant < self.probabilityMutant:
                                        a.status = Status.Mutant
                                        
                                self.schedule.add(a)
                                new_molecule_id += 1

        def step(self):
                self.status()
                self.move()
                self.age_increase()
                self.replicate()

model = Cell(30, 20, 20 ,0.01)
for index in range(20):
        model.step()
agent_state = model.datacollector.get_agent_vars_dataframe()
                

                


                

        
