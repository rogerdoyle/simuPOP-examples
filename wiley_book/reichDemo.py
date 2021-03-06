import math
class demoModel:
    def __init__(self, model, N0=1000, N1=100000, G0=500, G1=500, m=1):
        '''Return a demographic function with population split and expansion
        model: 'linear' or 'exponential'
        N0:   Initial population size.
        N1:   Ending population size.
        G0:   Length of burn-in stage.
        G1:   Length of population expansion stage.
        m:    Split population into m subpopulations before population expansion.
        '''
        self.N0, self.N1, self.G0, self.G1, self.m = N0, N1, G0, G1, m
        self.model = model
    
    def __call__(self, gen, pop=None):
        if self.model == 'instant':
            return self.ins_expansion(gen, pop)
        else:
            return self.exp_expansion(gen, pop)
    
    def ins_expansion(self, gen, pop):
        if gen < self.G0:
            return self.N0
        elif self.m > 1:
            if gen == self.G0:  # split population
                # avoid floating point problem
                pop.splitSubPop(0, [1./self.m] * self.m)
            return [self.N1//self.m]*(self.m-1)+[self.N1-self.N1//self.m*(self.m-1)] 
        else:
            return self.N1
    
    def exp_expansion(self, gen, pop):
        rate = (math.log(self.N1) - math.log(self.N0))/self.G1
        if gen < self.G0:
            return self.N0
        elif self.m > 1 and gen == self.G0:  # split population
            pop.splitSubPop(0, [1. / self.m] * self.m)
        if gen == self.G0 + self.G1 - 1:
            N = self.N1
        else:
            N = int(self.N0 * math.exp((gen - self.G0) * rate))
        if self.m > 1:
            return [N // self.m] * (self.m-1) + [N - N // self.m * (self.m-1)]
        else:
            return N

