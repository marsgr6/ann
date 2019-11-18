import random as r
from scipy.stats import binom
import warnings
import numpy as np

class TopologyException(Exception):
    def __init__(self, msg):
        super().__init__(msg)

class Topology:
    '''
    It describes the network topology. 

    Attributes
    ----------
    type : str
        Topology type: random, ring, or small-world. By default: random.

    omega : float
        Use together with small-world. By default omega=0.0.
        

    '''
    def __init__(self, type_="random", omega=1.0):
        self.type_=type_
        self.omega=omega

    @property 
    def type_(self):
        '''
        Get topology type
        '''
        return self._type_

    @type_.setter
    def type_(self, value):
        '''
        Set topology type
        '''
        topologies = ["random", "ring", "small-world"]
        if (type(value) != str) or (value not in topologies):
            raise TopologyException('The network topology should be a string from: [' +
                            ', '.join(topologies) + '].')
        self._type_ = value

    @property 
    def omega(self):
        '''
        Get omega
        '''
        return self._omega

    @omega.setter
    def omega(self, value):    
        '''
        Set omega
        '''
        if type(value) != float:
            raise ValueError('The type of the omega parameter '
                            'does not match the expected type (%s).'
                            % 'float')
        self._omega = value

    def make_topology(self, degree, neurons, _connectivity = "directed"):
        """
        Parameters
        ----------
        degree : 

        neurons : int
            Number of neurons
        _conectivity: str


        Returns
        -------
        A list of edges for each node

        """
        self.connectivity = _connectivity
        if self.type_ == "random":
            # Generates random network
            edges = self.random_network(degree, neurons, _connectivity)
        elif self.type_ == "ring":
            self.omega = 0.0
            if _connectivity == "directed" and self.type_ == "ring": 
                warnings.warn("Warning!. This implementation does not support the type of non-directed connectivity. By default it will create an undirected ring network.", Warning)
            # Generates ring network
            edges = self.ring_network(degree, neurons)        
        elif self.type_ == "small-world":
            self.omega = omega
            if _connectivity == "directed" and self.type_ == "small-world": 
                warnings.warn("Warning!. This implementation does not support the type of non-directed connectivity. By default it will create an undirected small world network.", Warning)
            if omega == 0.0:
                # Starts in a ring network
                edges = self.ring_network(degree, neurons)

            elif omega > 0.0 and omega < 1.0:
                edges = self.ring_network(degree, neurons)
                # Rewires ring network if w>0
                self.edges =self.rewireNet(self.omega, self.edges, self.neurons)

            elif omega == 1.0:
                edges = self.random_network(degree, neurons, _connectivity)
        return edges

    def generateEdge(self, ni, l, Nodes):
        """
        Generate new Edge for small-world network.

        Parameters
        ----------
        ni : type
            Description
        l : type
            Description

        Returns
        -------
        Description
        
        """
        while True:
            nE = r.randint(0, Nodes - 1)
            if nE != ni and nE not in l:
                return nE

    def rewireNet(self, p, C, Nodes):
        """
        Generates small-world network from
        regular ring, rewiring with p probability

        Parameters
        ----------
        p : float
            Rewiring probability
        C : list
            Network array 

        Returns
        -------
        
        """
        for ni in range(Nodes):
            for ci in C[ni]:
                if r.random() < p:
                    C[ni].remove(ci)  # Removing edge
                    C[ci].remove(ni)

                    nE = self.generateEdge(ni, C[ni], Nodes)  # Generates new Edge

                    C[ni].append(nE)  # Adding new edge
                    C[nE].append(ni)
        return C

    def random_directed_adjacency(self, N, p):
        """

        Parameters
        ----------

        Returns
        -------

        """
        for i in range(N):
            adjacency = binom.rvs(1, p, size=N)
            adjacency[i] = 0
            return np.where(adjacency)[0]

    def random_network(self, degree, nodes, _type="undirected"): 
        """
        Generates undirected Erdos-Renyi random network with
        edge creation probability p=Degree/Neurons=d/Nodes

        Parameters
        ----------
        degree : float
            Probability of the degree distribution 
        nodes : int
            Number of nodes 
        _type = str
            Conectivity

        Returns
        -------
        Array with an undirected Erdos-Renyi random network. 

        """
        p = float(degree) / nodes
    
        if _type == "undirected":
            C = [[] for x in range(nodes)]
            for i in range(nodes):
                for j in range(i+1, nodes):
                    if r.random() < p:
                        C[i].append(j)
                        C[j].append(i)
            
            return C
        
        if _type == "directed":
            return [self.random_directed_adjacency(nodes, p) for i in range(nodes)]
  

    def leftKneighbors(self, ni, d, Nodes):
        """

        Generates d/2 left neighbors

        Parameters
        ----------
        ni : int
            Description
        d  : float
            Description
        Nodes : int
            Description

        Returns
        -------
            A list of left network neightbors
        
        """
        lK = []
        for i in range(d // 2):
            ki = ni-i-1
            if ki < 0:
                lK.append(ki+Nodes)
            else:
                lK.append(ki)
            
        return lK


    def rightKneighbors(self, ni, d, Nodes):
        """

        Generates d/2 right neighbors

        Parameters
        ----------
        ni : int
            Description
        d  : float
            Description
        Nodes : int
            Description

        Returns
        -------
            A list of right network neightbors
        
        """
        rK = []
        for i in range(d // 2):
            ki = ni+i+1
            if ki < Nodes:
                rK.append(ki)
            else:
                rK.append(ki-Nodes)
            
        return rK

    def ring_network(self, degree, nodes):
        """
        Generates regular ring network with d/2 neighbors at each side
        
        Parameters
        ----------
        degree  : float
            Description
        nodes : int
            Description

        Returns
        -------
            A list with a ring network
        
        """
        return [self.leftKneighbors(ni, degree, nodes) + self.rightKneighbors(ni, degree, nodes) for ni in range(nodes)]