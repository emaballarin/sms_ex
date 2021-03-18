# Simulate SIR model

# Consider a system of differential equations $\dot{x} = f(x,t)$. If time is discretized, with a time interval $\Delta t$, the evolution in time can easily be computed as:
# $$
# x_{i+1} = x_i + f(x_i, t_i)\cdot \Delta t,
# $$
# where $t_i = t_0+i\cdot\Delta t$ and $x_i = x(t_i)$.

# Consider now the simple *SIR epidemiological model*, where $s$ is the number of susceptible population, $i$ is the number of infected and $r$ is the number of recovered/removed population. The population size $n$ is constant in time, $\forall t\ s(t)+i(t)+r(t)=n$. The variable $x(t) = (s(t), i(t), r(t))$.

# The system of ODE is
# $$
# \dot{x} =
# \begin{cases}
# -\beta\cdot s(t)\cdot i(t)\\
# \beta\cdot s(t)\cdot i(t)-\alpha\cdot i(t)\\
# \alpha\cdot i(t).\\
# \end{cases}
# $$

# Define a function that given an initial state $x(0)$ and a time $T$, simulates the evolution of the population until a time $T$ starting from $x(0)$.


import matplotlib.pyplot as plt
import numpy as np


def sir_ode(x, beta=0.1, alpha=0.01):

    return


# Plot the trajectories of the states of the system for different initial configurations.


# At each time step add some Gaussian noise to the state of the population. Multiply it by a constant $\sigma$ to govern the magnotude of the noise.

# use np.random.randn()
