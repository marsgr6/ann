import matplotlib.pyplot as plt
from ann import attractor_network as AA
from ann import pattern as p
import numpy as np

N = 1024 # Number of neurons
K = 4 # Degree
a = 0.5 # Sparseness
T = 800 # Steps


# Create network with same parameters as above
activity_net = AA.SimpleAttractor(N, K, a, "ring")
print(activity_net)
# Make topology (random by default)
activity_net.generate_topology()
# Make weights (random -1, 1)
activity_net.make_weights()
    # Initializa network state (binary random with a activity)
activity_net.make_initialization()
print(activity_net)v

t1 = datetime.datetime.now()
# Update network T steps (returns activity for each step)
activity = activity_net.update_steps(T)
t2 = datetime.datetime.now()
print(t2-t1)

# Plot activity
plt.figure(figsize=(20,6))
plt.plot(activity)
plt.show()