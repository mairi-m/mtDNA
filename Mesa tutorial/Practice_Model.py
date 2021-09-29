# Title: Practice Model
# Author: Mairi MacIain

from mesa import Agent, Model
from mesa.time import RandomActivation
import matplotlib.pyplot as plt
import random

# Create agent class
class TheAgent(Agent):
    def __init__(self, unique_id, model):
        self.unique_id = random.randint(1,100)
        print(unique_id)

    def step(self):
        return


class TheModel(Model):
    def __init__(self, N):
        self.num_agents = N
        self.schedule = RandomActivation(self)

        for index in range(self.num_agents):
            a = TheAgent(index, self)
            self.schedule.add(a)


    def step(self):
        self.schedule.step()

model = TheModel(5)
for index in range(5):
    model.step()

        
        
    
    











        
        
