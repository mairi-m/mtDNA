# Title: Mesa intro
# Author: Mairi MacIain
# Date: 25 June 2021

from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.space import MultiGrid
import matplotlib.pyplot as plt
import numpy as np


class MoneyAgent(Agent):
        # An agent with fixed initial wealth
        def __init__(self, unique_id, model):
                super().__init__(unique_id, model)
                self.wealth = 1


        def move(self):
                possible_steps = self.model.grid.get_neighborhood(
                        self.pos,
                        moore = True,
                        include_center = False)

                new_position = self.random.choice(possible_steps)
                self.model.grid.move_agent(self, new_position)


        def give_money(self):
                cellmates = self.model.grid.get_cell_list_contents([self.pos])
                if len(cellmates) > 1:
                        other = self.random.choice(cellmates)
                        other.wealth += 1
                        self.wealth -= 1


        def step(self):
                self.move()
                if self.wealth > 0:
                        self.give_money()
                


class MoneyModel(Model):
        # A model with some number of agents
        def __init__(self, N, width, height):
                self.num_agents = N
                self.grid = MultiGrid(width, height, True)
                self.schedule = RandomActivation(self)

                # Create agents
                for index in range(self.num_agents):
                        a = MoneyAgent(index, self)
                        self.schedule.add(a)

                        # Add the agent to a random grid cell
                        x = self.random.randrange(self.grid.width)
                        y = self.random.randrange(self.grid.height)
                        self.grid.place_agent(a, (x, y))

        def step(self):
                # Advance the model by one step
                self.schedule.step()

                
model = MoneyModel(50, 10, 10)
for index in range(20):
        model.step()


agent_counts = np.zeros((model.grid.width, model.grid.height))
for cell in model.grid.coord_iter():
        cell_content, x, y = cell
        agent_count = len(cell_content)
        agent_counts[x][y] = agent_count
plt.imshow(agent_counts, interpolation = 'nearest')
plt.colorbar()

plt.show()



