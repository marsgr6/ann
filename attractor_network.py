import numpy as np

from .topology import Topology
from .pattern import Pattern
from .network import Network
from scipy.stats import binom

class Attractor:
    '''
    Description : Network attractor

    Attributes
    ----------
    neurons : int
        Number of neurons.
    degree : int
        Degree distribution.
    sparseness : type
        It is the proportion of neurons active to represent any one stimulus in the network. 
    topology : type
        It refers to the network' topology ( "random", "ring", "small-world" ).
    omega : type
        Random shortcut probability
    activity : float
        Neural activity
    theshold : float
         It keeps  network activity stay at expected values ​​close to the sparseness parameter

    '''
    def __init__(self, nodes, degree, sparseness, topology='undefined', omega=None):
        self.neurons = nodes
        self.degree = degree
        self.sparseness = sparseness
        self.threshold = (1-2*self.sparseness)/(
            2*np.sqrt(self.sparseness*(1-self.sparseness)))
        self.topology = topology
        self.omega = omega
        self.activity = 0.0
        self.network = Network() 

    @property 
    def neurons(self):
        '''
        Get neurons
        '''
        return self._neurons

    @neurons.setter
    def neurons(self, value):
        '''
        Set neurons
        '''
        if type(value) != int:
            raise ValueError('The type of the nodes parameter '
                        'does not match the expected type (%s).'
                        % 'int')
        self._neurons = value

    @property
    def degree(self):
        '''
        Get degree
        '''
        return self._degree

    @degree.setter
    def degree(self, value):
        '''
        Set degree
        '''
        if type(value) != int:
            raise ValueError('The type of the degree parameter '
                        'does not match the expected type (%s).'
                        % 'int')
        self._degree = value

    @property
    def sparseness(self):
        '''
        Get sparseness
        '''
        return self._sparseness
    
    @sparseness.setter
    def sparseness(self, value):
        '''
        Set sparseness
        '''
        if type(value) != float:
            raise ValueError('The type of the sparseness parameter '
                        'does not match the expected type (%s).'
                        % 'float')
        self._sparseness = value

    def random_weights(self, C):
        """
        Generate random weights for the network

        Parameters
        ----------
        C : list
            List of network edges 

        Returns
        -------
        A list with the updated network  weights

        """
        return [binom.rvs(1, 0.5, size=len(c)) * 2 - 1 for c in C]

    def zero_weights(self, C):
        """
        Set the network weiths to 0

        Parameters
        ----------
        C : lsit 
            List of network edges

        Returns
        -------
        A list with the updated network  weights

        """
        return [np.zeros(len(c), dtype='float') for c in C]

    def make_weights(self, _type="random"):
        """
        Random weights from {-1, 1}

        Parameters
        ----------
        _type: str
            Weight generation type

        """
        if _type == "random": 
            weights = self.random_weights(self.network.get_edge_list())
        if _type == "zero":
         weights = self.zero_weights(self.network.get_edge_list())
        self.network.update_edge_weights(weights)

    def sign_function(self, s_value):
        """
        
        Parameters
        ----------
        s_value : float 
            Description 

        Returns
        -------
        Argument value [-1, 1]
                
        """
        if s_value >= 0: 
            return 1
        return -1

    def heavyside_function(self, s_value):
        """

        Unit step function where the arguments can be zero for negative and one for positive.

        Parameters
        ----------
        s_value : float 
            Description 

        Returns
        -------
        Argument value [0, 1]
        
        """
        if s_value >= 0:
            return 1
        return 0
    

    def synchronous_update(self, _type = "heavyside"):
        """

        Synchronous update of neurons states

        Parameters
        ----------

        Returns
        -------

        """
        state_new = np.zeros(self.neurons, dtype='int')
        for idx in self.network.nodes:
            node_id = self.network.nodes[idx].key
            f_value = np.sum( self.network.get_node_neighbors_weights(node_id) 
                * self.network.get_node_neighbors_states(node_id))
            if _type == "heavyside": 
                state_new[node_id] = self.heavyside_function(f_value)
            elif _type == "sign": 
                state_new[node_id] = self.sign_function(f_value)
        # Neural activity is updated
        self.network.update_states(state_new)
        self.activity = np.mean(state_new)

    def make_initialization(self, type_="random"):
        """

        Binary initial state with given sparseness

        Parameters
        ----------
        type_ : str
            Network type
        
        """
        if type_ == "random":
            state_new = binom.rvs(1, self.sparseness, size=self.neurons)
            self.network.update_states(state_new)
        self.activity = np.mean(self.network.get_states())

    def generate_topology(self, _type="random", omega = 1.0):
        """
        It generates a topology based on a omega value.

        Return
        -------
        
        """
        if self.topology == 'undefined':
            t = Topology()
        elif self.topology !='undefined' and self.omega==None: 
            t = Topology(self.topology)
        else:
            t = Topology(self.topology, self.omega)
        self.topology = t.type_
        self.omega = t.omega
        self.network.generate_from_edge_list(t.make_topology(self.degree, self.neurons))

    def __str__(self): 
        """

        Print the global network parameters

        Returns
        -------
        String with the global network parameters
        
        """
        return "Network with: nodes={}, degree={}, sparseness={}, topology={}, omega={},  theta={}".format(
            self.neurons, self.degree, self.sparseness, self.topology, str(self.omega), str(self.threshold))       

class SimpleAttractor(Attractor):

    def update_steps(self, steps, metric="activity"):
        """
        

        Parameters
        ----------
        steps : int
            Number of steps
            
        metric: str
            Metric activity type
    
        Returns
        -------
        
        """
        if metric == "activity":
            a = []
            for t in range(steps):
                self.synchronous_update() 
                a += [self.activity]
            return a

        if metric == "overlap":
            pass

class ClassicAttractor(Attractor):
    """
    Classic attractor with polar states {-1,1}

    See Also
    --------
    https://neuronaldynamics.epfl.ch/online/Ch17.S2.html

    """

    def update_weights(self):
        """
        Hebbian learning

        See Also
        --------
        https://en.wikibooks.org/wiki/Artificial_Neural_Networks/Hebbian_Learning

        """
        
        """
        weights=[]
        for idx in sorted(self.network.nodes):
            node_id = self.network.nodes[idx].key
            weights +=  [self.network.get_node_neighbors_states(node_id) * self.network.get_node_state(node_id)] #/ self.neurons
        self.network.update_edge_weights(weights)
        """
        edges = self.network.get_edge_list()
        weights = self.network.get_edge_weights()
        for ci, adj_ci in enumerate(edges):
            weights[ci] += self.network.get_node_neighbors_states(ci) * self.network.get_node_state(ci) #/ self.neurons
        self.network.update_edge_weights(weights)

        

    def update_steps(self, steps, pattern):
        """
        
        Parameters
        ----------
        steps : int
        pattern : int 

        Returns
        -------
        
        """       
        a = [np.mean(self.network.get_states())]
        m = [np.sum(self.network.get_states() * pattern) / self.neurons]
        for t in range(steps):
            self.synchronous_update("sign")
            a += [self.activity]
            m += [np.sum(self.network.get_states()* pattern) / self.neurons]
            if m[-1] == m[-2]:
                break
        return a, m


class SparseAttractor(Attractor):      


    def update_weights(self):
        """
        Normalize pattern before Hebb learning

        Parameters
        ----------

        
        """     
        a = np.mean(self.network.get_states())
        A = a * (1-a)
        normalized_state = (self.network.get_states() - a) / np.sqrt(A)
        normalized_state_list=list(normalized_state)
        edges = self.network.get_edge_list()
        weights = self.network.get_edge_weights()
        for ci, adj_ci in enumerate(edges):
            weights[ci] += np.array([normalized_state_list[x.key] for x in adj_ci]) * normalized_state_list[ci]#self.network.get_node_neighbors_states(ci) * self.network.get_node_state(ci) #normalized_state[adj_ci] * normalized_state[ci] #/ self.neurons
        self.network.update_edge_weights(weights)      

    def update_steps(self, steps, pattern):
        """

        Parameters
        ----------

        Returns
        -------
        
        """        

        a_net = np.mean(self.network.get_states())
        A_net = a_net * (1-a_net)
        normalized_state = (self.network.get_states() - a_net) / np.sqrt(A_net)
        
        a_pat = np.mean(pattern)
        A_pat = a_pat * (1-a_pat)
        normalized_pattern = (pattern - a_pat) / np.sqrt(A_pat)
        
        activity = [a_net]
        overlap = [np.sum(normalized_state * normalized_pattern) / self.neurons]
        
        for t in range(steps):
            self.synchronous_update()
            
            a_net = np.mean(self.network.get_states())
            A_net = a_net * (1-a_net)
            normalized_state = (self.network.get_states() - a_net) / np.sqrt(A_net)
        
            activity += [self.activity]
            overlap += [np.sum(normalized_state * normalized_pattern) / self.neurons]
            
            if overlap[-1] == overlap[-2]:
                break
        return activity, overlap, t
    
    
    def synchronous_update(self):
        """
        Synchronous update of neurons states

        Parameters
        ----------
        
        """        
        state_new = np.zeros(self.neurons, dtype='int')
        
        for i in self.network.nodes:
        #for i in range(self.neurons):
            neighborhood_i = self.network.get_node_neighbors_states(i)#self.state[self.edges[i]]
            a = np.mean(neighborhood_i)
            A = a * (1-a)
            neighborhood_i_normalized = (neighborhood_i - a) / np.sqrt(A)
            f_value =  np.sum(self.network.get_node_neighbors_weights(i) * neighborhood_i_normalized) /self.degree#np.sum(self.weights[i] * neighborhood_i_normalized) / self.degree
            f_value = f_value - self.threshold
            state_new[i] = self.heavyside_function(f_value)
        
        self.network.update_states(state_new)
        
        self.activity = np.mean(self.network.get_states())
        


class EnsembleAttractor():
    '''
    Ensemble Attractor

    Attributes
    ----------
    neurons : int
        Number of neurons.
    modules
    connectivity
    sparseness
    topology
    
    degree : int
        Degree distribution.
    sparseness : type
        It is the proportion of neurons active to represent any one stimulus in the network. 
    topology : type
        It refers to the network' topology ( "random", "ring", "small-world" ).
    omega : type
        Random shortcut probability
    activity : float
        Neural activity
        
    '''

    def __init__(self, neurons, modules, connectivity, sparseness, topology, threshold = None):
        self.anns = []
        self.K = np.sum(connectivity)
        self.create_modules(neurons, modules, connectivity, sparseness, topology, threshold)



    def create_modules(self, neurons, modules, connectivity, sparseness, topology, threshold): 
        """
        Description

        Parameters
        ----------
        neurons : int
            Description
        modules : 
            Description 
        connectivity :
            Description
        sparseness :
            Description
        topology :
            Description
        threshold : 
            Description

        """
        for i in range(modules):
            self.anns += [SparseAttractor(neurons, connectivity[i], sparseness[i])]
            self.anns[i].generate_topology(topology[i][0], omega=topology[i][1])
            self.anns[i].make_weights("zero")
            self.anns[i].make_initialization("zero")
            if threshold != None:
                self.anns[i].threshold = threshold[i]
            print(self.anns[i])        


    def train_and_test_offline(self, train_set, test_set):
        """

        Parameters
        ----------
        train_set : array
            Description
        test_set : array
            Description

        """
        M = []
        A = []
        T = []
        
        # Train
        for si, set_i in enumerate(train_set):
            for pattern_i in set_i:
                self.anns[si].network.update_states(pattern_i) 
                self.anns[si].update_weights()
        
        # Test
        for si, set_i in enumerate(test_set):
            for pattern_i in set_i:
                self.anns[si].network.update_states(pattern_i) 
                a, m, ts = self.anns[si].update_steps(100, pattern_i)
                M += [m[-1]]
                A += [a[-1]]
                T += [ts]
                
        return np.array(M)


    def train_and_test_online(self, train_set, test_set, mu, theta):
        """

        Parameters
        ----------
        train_set : array
            Description
        test_set : array
            Description
        mu : float
            Description
        theta : 
            Description

        """        
        mean_O = []
        R = []
        alpha_R = []
        modules = len(train_set)
        
        # Incrementaly learn and test
        for p_i in range(mu):
            # Train
            for si, set_i in enumerate(train_set): 
                self.anns[si].network.update_states(set_i[p_i]) 
                self.anns[si].update_weights()
                
            # Test
            M = []
            A = []
            T = []
            for si, set_i in enumerate(test_set):
                for tp_i in range(p_i + 1):
                    self.anns[si].network.update_states(set_i[p_i])
                    a, m, ts = self.anns[si].update_steps(100, set_i[tp_i])
                    M += [m[-1]]
                    A += [a[-1]]
                    T += [ts]
                    
            M = np.array(M)
            mean_O += [np.mean(M)]
            P_r = len(np.where(M > theta)[0])
            P_l = modules * (p_i + 1)
            R += [P_r / P_l]
            alpha_R += [P_r / self.K]
            
            
        return np.array(mean_O), np.array(R), np.array(alpha_R)
                    
