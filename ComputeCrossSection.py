import numpy as np


class ComputeCrossSection(object):
    """Documentation for ComputeCrossSection

    """

    def __init__(self, ray, pol):
        super(ComputeCrossSection, self).__init__()
        self.ray = ray
        self.pol = pol
        self.cd = np.block([np.matrix([0 for i in range(np.shape(pol)[1])]),
                            np.matrix([0])])
        self.cb = np.matrix([1, 1, 1])
        self.A = np.block([[pol, -self.ray.T],
                           [np.matrix([1 for i in range(np.shape(pol)[1])]),
                            np.matrix([0])]])
        self.B = np.matrix([[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]])
        self.b = np.matrix([[self.ray0[0, 0], self.ray0[0, 1], 1]])
        self.Binv = np.linalg.inv(self.B)
        for i in range(np.shape(self.b)[1]-1):
            if self.b[0, i] < 0:
                self.b[0, i] *= -1
                self.A[i] *= -1
        self.D = self.A
