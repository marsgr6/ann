# Periodic Neural Activity

Periodic Neural Activity Induced by Network Complexity

## Overview


## Requirements

Python 3.5+

## Quick Start

```
import matplotlib.pyplot as plt
from an_library import attractor_network as AA
from an_library import pattern as p

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
