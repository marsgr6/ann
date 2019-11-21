import matplotlib.pyplot as plt
from ann import attractor_network as AA
from ann import pattern as p
import numpy as np

# Network parameters
N = 1000
K = 100
a = 0.5

# Create network, topology, weights and initial state
classic_net = AA.ClassicAttractor(N, K, a)
print(classic_net)
classic_net.generate_topology()
classic_net.make_weights("zero")
classic_net.make_initialization("zero")
print(classic_net)

# Learn P patters
P = 65
M = []  # Overlap vector for each pattern
A = []  # Activity vector for each pattern
for i in range(P):
    print("\r%d" % (i+1), end='')
    pattern = p.Pattern.generate_random_pattern(classic_net.neurons, _type="polar")  # Generate pattern
    classic_net.network.update_states(pattern)  # Make network state equal to the pattern
    classic_net.update_weights()  # Update weights (Hebbian Learning)
    a, m = classic_net.update_steps(100, pattern)  # Update network 100 times
    M += [m[-1]]  # Save last value of m for the pattern (after 100 times)
    A += [a[-1]]  # Sava last value of a for the pattern (after 100 times

plt.figure(figsize=(14,7))

plt.figure(1)
plt.subplot(121)
plt.axis([0, P/K, 0, 1])
alpha = np.array(range(1,P+1)) / K
plt.plot(alpha, np.abs(M), ':o')
plt.subplot(122)
plt.plot(A)