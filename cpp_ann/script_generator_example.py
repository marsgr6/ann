"""
  This code generates a unix bash script  
  for running the C++ generated program 
  with the selected set of parameters for 
  p: patterns (integer)
  K: connectivity (integer)
  nns: number of modules (list of integers)
  topos: topologies (list from ['r', 'c', 'x'])
    r: Ring, c: Cross, 'x': X topologies
  ws: omega values (list of floats)

  Example:
    python script_generator_example.py > run_example.sh

  Then make run_example.sh executable: 
    chmod +x run_example.sh

  Execute:
    ./run_example.sh

  A list of files for each execution of the C++ program will be found in 
  the . directory.
"""
p = 150  # patterns
K = 240

nns = [1, 10, 15]  # 10 and 15 ensemble modules
topos = ['r', 'c']  # r: ring topology
ws = [0.0, 0.5, 1]  #

advance = 0

runs = sum([p // nn * len(topos) * len(ws) for nn in nns])

path = "patterns/"

# Iterates for each list above
for nn in nns:
    for topo in topos:
        for w in ws:
            for i in range(p // nn):
                advance += 1
                percentage = advance / runs * 100
                print(r"echo -e '\r" + "%0.0f" % percentage + "%'")
                print("./sparsenet 89420 " +
                      str(K // nn) + " " +
                      str(w) +
                      " 0.2258 1 r 0.656 0.7 0.0 100 1 " +
                      str((i + 1) * nn) +
                      " 6 100 263 340 " + path + " " + path + " " + 
                      topo + " " + str(i + 1) + " " + str(nn))
