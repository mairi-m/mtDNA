# Title: Infection model

# Imports
import time
import numpy as np
import pandas as pd
import pylab as plt
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector

class InfectionModel(Model):
    # A model for infection spread

    def __init__(self, N = 10, width = 10,height = 10, ptrans = 0.5, death_rate = 0.02, recovery_days = 21, recovery_sd = 7):
        
        self.num_agents = N
        self.recovery_days = recovery_days
        self.recovery_sd = recovery_sd
        self.ptrans = ptrans
        self.death_rate = death_rate
        self.grid = MultiGrid(width, height, True)
        self.schedule = RandomActivation(self)
        self.running = True
        self.dead_agents = []

        # Create agents

        for index in range(self.num_agents):
            a = person(index, self)
            self.schedule.add(a)

            # Add the agent to a random grid cell

            x = self.random.randrange(self.grid.width)
            y = self.random.randrange(self.grid.height)
            self.grid.place_agent(a, (x,y))

            # Make some agents infected at the start
            infected = np.random.choice([0,1], p=[0.98,0.02])
            if infected == 1:
                a.state = State.INFECTED
                a.recovery_time = self.get_recovery_time()

        self.datacollector = DataCollector(
            agent_reporters = {"State": "state"})

    def get_recovery_time(self):
        return int(self.random.normalvariate(self.recovery_days.recovery_sd))


    def step(self):
        self.datacollector.collect(self)
        self.schedule.step()

class State(enum.IntEnum):
    SUSCEPTIBLE = 0
    INFECTED = 1
    REMOVED = 2

class person(Agent):
    # An agent in an epidemic model

    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.age = self.random.normalvariate(20,40)
        self.state = State.SUSCEPTIBLE
        self.infection_time = 0


    def move(self):
        # Move the agent

        possible_steps = self.model.grid.get_neighborhood(
            self.pos,
            moore = True,
            include_center = False)
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)

    def status(self):
    # Check infection status

        if self.state == State.INFECTED:
            drate = self.model.death_rate
            alive = np.random.choice([0,1], p=[drate, 1-drate])
            if alive == 0:
                self.model.schedule.remove(self)
            t = self.model.schedule.time-self.infection_time
            if t >= self.recovery_time:
                self.state = State.REMOVED

    def contact(self):
        # FInd close contacts and infect

        cellmates = self.model.grid,get_cell_list_contents([self.pos])
        if len(cellmates) > 1:
            for other in cellmates:
                if self.random.random() > model.ptrans:
                    continue
                if self.state is State.INFECTED and other.state is State.SUSCEPITBLE:
                    other.state = State.INFECTED
                    other.infection_time = self.model.schedule.time
                    other.recovery_time = self.model.get_recovery_time()

    def step(self):
        self.status()
        self.move()
        self.contact()



# RUN THE MODEL

model = InfectionModel(pop, 20, 20, ptrans = 0.5)
for i in range(steps):
    model.step()
agent_state = model.datacollector.get_agent_vars_dataframe()


def get_column_data(model):
    # Get states count at each step
    agent_state = model.datacollector.get_agent_vars_dataframe()
    X = pd.pivot_table(agent_state.reset_index(), index = 'Step', columns = 'State', aggfunc = np.size, fill_value = 0)
    labels = ['Susceptible', 'Infected', 'Removed']
    X.columns = labels[:len(X.columns)]
    return X








        
