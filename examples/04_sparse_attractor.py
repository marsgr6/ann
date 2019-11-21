import matplotlib.pyplot as plt
from ann import attractor_network as AA
from ann import pattern as p
import numpy as np

N = 1000
K = 100
a = 0.5

sparse_net = AA.SparseAttractor(N, K, a)
print(sparse_net)
sparse_net.generate_topology()
sparse_net.make_weights("zero")
sparse_net.make_initialization("zero")
print(sparse_net)



P = 60
M = []
A = []
T = []
for i in range(P):
    print("\r%d" % (i+1), end='')
    pattern = p.Pattern.generate_random_pattern(sparse_net.neurons)
    sparse_net.network.update_states(pattern) 
    sparse_net.update_weights()
    a, m, ts = sparse_net.update_steps(100, pattern)
    M += [m[-1]]
    A += [a[-1]]
    T += [ts]
    
plt.figure(figsize=(16,6))

plt.figure(1)
plt.subplot(131)
plt.axis([0, 1, 0, 1])
alpha = np.array(range(1,P+1)) / K
plt.plot(alpha, M, ':o')
plt.subplot(132)
plt.plot(A)
plt.subplot(133)
plt.plot(T)