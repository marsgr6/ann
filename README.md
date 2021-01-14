# Periodic Neural Activity
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/Naereen/StrapDown.js/blob/master/LICENSE)


## Overview


## Requirements

Python 3.5+

## Quick Start

Code below replicates the base result of 
the paper Periodic Neural Activity Induced by Network Complexity, https://journals.aps.org/pre/abstract/10.1103/PhysRevE.74.017102. 

```
import matplotlib.pyplot as plt
from ann import attractor_network as AA
from ann import pattern as p

N = 1024 # Number of neurons
K = 4 # Degree
a = 0.5 # Sparseness
T = 800 # Steps

# Create network
activity_net = AA.SimpleAttractor(N, K, a)
# Make topology (random by default)
activity_net.generate_topology()
# Make weights (random -1, 1)
activity_net.make_weights()
# Initializa network state (binary random with a activity)
activity_net.make_initialization()

# Update network T steps (returns activity for each step)
activity = activity_net.update_steps(T)

# Plot activity
plt.figure(figsize=(20,6))
plt.plot(activity)
plt.show()
```
