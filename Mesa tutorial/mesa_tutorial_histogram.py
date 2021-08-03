# Title: Mesa intro
# Author: Mairi MacIain
# Date: 25 June 2021

from mesa import Agent, Model
from mesa.time import RandomActivation
import matplotlib.pyplot as plt

class MoneyAgent(Agent):
        # An agent with fixed initial wealth
        def __init__(self, unique_id, model):
                super().__init__(unique_id, model)
                self.wealth = 1

        def step(self):
                if self.wealth == 0:
                        return
                other_agent = self.random.choice(self.model.schedule.agents)
                other_agent.wealth += 1
                self.wealth -= 1


class MoneyModel(Model):
        # A model with some number of agents
        def __init__(self, N):
                self.num_agents = N
                self.schedule = RandomActivation(self)

                # Create agents
                for index in range(self.num_agents):
                        a = MoneyAgent(index, self)
                        self.schedule.add(a)

        def step(self):
                # Advance the model by one step
                self.schedule.step()

all_wealth = []
for j in range(100):
    
    model = MoneyModel(50)
    for index in range(10):
        model.step()

    for agent in model.schedule.agents:
        all_wealth.append(agent.wealth)

plt.hist(all_wealth, bins=range(max(all_wealth)+1))




plt.show()

