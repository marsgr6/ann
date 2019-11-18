from scipy.stats import binom

class Pattern:
    '''

    '''
    def generate_random_pattern(p_size, activity = 0.5, _type = "binary"):
        """

        Parameters
        ----------
        p_size: int
            Number of trials.
        activity : int 
            It denotes the probability of success in each trial.
        _type: str
            It refers to the pattern type( "binary", "polar" ). 

        Returns
        -------
        It returns a vector with n trials 
        
        """    
        pattern = binom.rvs(1, activity, size=p_size)
        
        if _type == "polar":
            pattern = pattern*2 - 1
        return pattern


    def generate_random_pattern_set(p_size, patterns, activity = 0.5, _type = "binary"):
        '''

        Parameters
        ----------
        p_size : int

        patterns : array

        activity : float

        _type: str


        Returns
        -------

        '''
        pattern_set = binom.rvs(1, activity, size=(patterns, p_size))
        
        if _type == "polar":
            pattern_set = pattern_set*2 - 1
        
        return pattern_set
