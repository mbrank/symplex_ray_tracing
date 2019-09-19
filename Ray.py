class Ray(object):
    """Ray class

    """
    def __init__(self, R0, R):
        super(Ray, self).__init__()
        self.R0 = R0
        self.R = R
        self.Lam = Lam


    def setR0(self, R0):
        self.R0 = R0
        
        
    def getR0(self):
        return self.R0


    def setR(self, R):
        self.R = R
        
        
    def getR(self):
        return self.R
    

    def setLam(self, Lam):
        self.Lam = Lam
        
        
    def getLam(self):
        return self.Lam
    
