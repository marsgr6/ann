import numpy as np

class Neightbor(object):
    '''
    Neighbor node. It's associate to a network node

    Attributes
    ----------
    key : int 
        Node identifier
    weight : float
        It refers to the edge weight

    '''

    def __init__(self, key, weight=0):
        self.key = key
        self.weight=weight

    def __repr__(self):
        '''
        Print the neightbor detail
        '''
        return '{}'.format(self.key)

class Neuron(object):
    '''
    Description 

    Attributes
    ----------
    key : int
        identifier of the neuron 
    state: int
        Firing (1=active) of non-firing(0=inactive)

    '''

    def __init__(self, key, state=0):
        self.key = key
        self.neighbors = {}
        self.state=state

    def add_neighbor(self, neighbor):
        """
        Add a neighbor to a network node

        Parameters
        ----------
        neighbor: object
            It is a neuron
        
        """    
        id_=len(self.neighbors)
        self.neighbors[id_] = Neightbor(neighbor)


    def get_connections(self):
        """
        Get the neuron conections 

        Returns
        -------
        A list with the neuron neighbors
            
        """    
        return self.neighbors.keys()

    def get_weight(self, idx):
        """

        Parameters
        ----------
        idx : int
            Neighbor identifier

        Returns
        -------
        It returns an edge weight

        """    
        return self.neighbors[idx].weight

    def __repr__(self):
        return '{}({}) -> {}'.format(
            self.key, self.state, [self.neighbors[x].key for x in self.neighbors])


class Network(object):
    '''
    Neuron network 

    Attributes
    ----------
    nodes : dic
        List of neurons

    '''
    
    def __init__(self):
        self.nodes = {}


    def get_vertices(self):
        """
        Get the ids of the network neurons

        Returns
        -------
        Neuron list
        
        """    
        return self.nodes.keys()


    def add_vertex(self, vertex):
        """
        Adds a neuron to the network

        Parameters
        ----------
        vertex : obj
            Neuron
        
        """    
        self.nodes[vertex.key] = vertex

    def get_vertex(self, key):
        """
        Get a neuron of the network

        Parameters
        ----------
        key : int
            Neuron identifier

        Returns
        -------
        A neuron

        """    
        try:
            return self.nodes[key]
        except KeyError:
            return None


    def add_edge(self, from_key, to_key):
        """

        Add an edge to the network

        Parameters
        ----------
        from_key : int
            Source identifier
        to_key : int
            Target identifier (Neighbor)
        
        """    
        if from_key not in self.nodes:
            self.add_vertex(Neuron(from_key))
        if to_key not in self.nodes:
            self.add_vertex(Neuron(to_key))
        self.nodes[from_key].add_neighbor(Neightbor(to_key))


    def generate_from_edge_list(self, edge_list): 
        """
        Generate a network from a list of list of edges neuron ids

        Parameters
        ----------
        edge_list : list
            A list of list with the ids of the neurons 
        """    
        [[self.add_edge(idx_a, idx_b) for idx_b in edge_list] for idx_a, edge_list in enumerate(edge_list)]


    def update_state(self, node_idx, new_node_state): 
        """
        Update a neuron state

        Parameters
        ----------
        node_idx : int
            Neuron identifier
        new_node_state: int
            Neuron state. Firing (1=active) of non-firing(0=inactive)
        """    
        self.nodes[node_idx].state=new_node_state

    def get_node_state(self, idx):
        """
        Get the state of a neuron

        Parameters
        ----------
        idx : int
            Neuron identifier

        Returns
        -------
        A neuron state

        """    
        return self.nodes[idx].state

    def update_states(self, state_list): 
        """
        Update states in the network

        Parameters
        ----------
        state_list : list
            A list of the states to be updated in the current network
        
        """    
        [self.update_state(node_idx, state_list[node_idx]) for node_idx in sorted(self.nodes)]

    def get_states(self): 
        """
        Get the states of the network

        Returns
        -------
        A list with the states of the network
        
        """    
        return [self.nodes[x].state for x in sorted(self.nodes)]

    def get_edge_list(self):
        """
        Get a list of edges for each neuron

        Returns
        -------
        A list of list with the ids of the neurons 

        """    
        return [[self.nodes[n_a].neighbors[n_b].key for n_b in self.nodes[n_a].neighbors ] for n_a in sorted(self.nodes)] 


    def get_edge_weights(self):
        """
        Get a list of list with the weights of each edge

        Returns
        -------
        A list of list with the edge' weights
        
        """    
        return [[self.nodes[idx_a].neighbors[idx_b].weight for idx_b in self.nodes[idx_a].neighbors] for idx_a in sorted(self.nodes)]


    def get_neighbor_id(self, from_key, to_key):
        for x in self.nodes[from_key].neighbors:
            neighbor=self.nodes[from_key].neighbors[x]
            if str(neighbor.key)==str(to_key):
                return x

    def update_edge_weight(self, from_key, to_key, weight_):
        """
        Update an edge weight

        Parameters
        ----------
        from_key : int
            Source identifier
        to_key : int
            Target identifier
        weight_ : float
            Edge weight
        """    
        self.nodes[from_key].neighbors[self.get_neighbor_id(from_key, to_key)].weight=weight_

    def update_edge_weight_with_position(self, from_key, position, weight_):
        """
        Update an edge weight

        Parameters
        ----------
        from_key : int
            Source identifier
        to_key : int
            Target identifier
        weight_ : float
            Edge weight
        """    
        self.nodes[from_key].neighbors[position].weight=weight_


    def update_edge_weights(self, edge_weight_list): 
        """
        Update the weights in the network

        Parameters
        ----------
        edge_weight_list : list
            A list of list with the edge' weights
        
        """    
        for i, current_edge_weigth in enumerate(edge_weight_list):
            for j, weight in enumerate(current_edge_weigth):
                self.update_edge_weight_with_position(i, j, weight)


    def get_node_neighbors_weights(self, idx):
        """
        Get the weights of a neuron' neighbors

        Parameters
        ----------
        idx : int 
            Neuron identifier

        Returns
        -------
        An array with the weights of the neuron's neighbors
        
        """    
        return np.array([self.nodes[idx].neighbors[x].weight for x in self.nodes[idx].neighbors])


    def get_node_neighbors_states(self, idx):
        """
        Get the ids of a neuron' neighbors

        Parameters
        ----------
        idx : int 
            Neuron identifier

        Returns
        -------
        An array with the ids of the neuron's neighbors
        
        """    
        return(np.array([self.nodes[int(str(self.nodes[idx].neighbors[x].key))].state for x in self.nodes[idx].get_connections()]))

    def __contains__(self, key):
        return key in self.nodes

    def __iter__(self):
        """
        Iterate between the network neurons

        Returns
        -------
        Nueron iterator
        
        """    
        return iter(self.nodes.values())


    def __repr__(self):
        """
        Print a list with the network neurons
        
        """    
        return "{}".format([self.nodes[x] for x in self.nodes])

    def print_network(self):
        """
        Print all the network detail
        node_id ( node status ) >> [neighbor_id : neighbor_weight, ...]
        
        """    
        for x in sorted(self.nodes):
            list_neightbors=[str(self.nodes[x].neighbors[n].key)+" : "+str(self.nodes[x].neighbors[n].weight) for n in self.nodes[x].neighbors]
            print("{}({}) >> {}".format(self.nodes[x].key, self.nodes[x].state, list_neightbors))



