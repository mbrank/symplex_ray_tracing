import numpy as np


class ComputeCrossSection(object):
    """Documentation for ComputeCrossSection

    """

    def __init__(self, ray, polynomial):
        super(ComputeCrossSection, self).__init__()
        self.ray = ray.getR()
        self.ray0 = ray.getR0()
        self.pol = polynomial.getListofPoints()

        self.cd = np.block([np.matrix([0 for i in
                                       range(np.shape(self.pol)[1])]),
                            np.matrix([0])])
        self.cb = np.matrix([1, 1, 1])
        self.A = np.block([[self.pol, -self.ray.T],
                           [np.matrix([1 for i in
                                       range(np.shape(self.pol)[1])]),
                            np.matrix([0])]])
        self.Binv = np.matrix([[1, 0, 0],
                               [0, 1, 0],
                               [0, 0, 1]])
        self.b = np.matrix([[self.ray0[0, 0], self.ray0[0, 1], 1]])
        self.Binv = np.linalg.inv(self.B)

        for i in range(np.shape(self.b)[1]-1):
            if self.b[0, i] < 0:
                self.b[0, i] *= -1
                self.A[i] *= -1
        self.D = self.A
        self.obj_func = ['z']
        self.real_vars = ['x'+str(i) for i in range(np.shape(self.pol)[1]+1)]
        self.lam = self.real_vars[-1]
        self.pivot_el_col = ['ax'+str(i) for i in range(np.size(self.b))]
        self.constraint = ['b']

    def checkIntersection(self):
        self.sigma = -self.cb*self.Binv*self.A
        self.xb = self.Binv*self.b.T
        self.z = self.cb*self.Binv*self.b.T
        self.z = self.solveMatrix(self.cb)
        if self.z[0, 0] != 0:
            print("None")
            return None
        else:
            self.cb = np.matrix([1, 1, 1])
            self.cd = np.block([np.matrix([0 for i in
                                           range(np.shape(self.pol)[1])]),
                                np.matrix([1])])
            self.solveMatrix(self.cb)
            res_dict = {}
            for i in range(len(self.pivot_el_col)-1):
                res_dict[self.pivot_el_col[i]] = np.asarray(
                    self.Binv*self.b.T)[i, :][0]

            lam = self.real_vars[-1]
            if lam in res_dict:
                ex1 = self.ray0[0, 0]+self.ray[0, 0]*res_dict[lam]
                ey1 = self.ray0[0, 1]+self.ray[0, 1]*res_dict[lam]
                import matplotlib.pyplot as plt
                plt.scatter(ex1, ey1)
                plt.text(ex1+0.1, ey1+0.1, "e", fontsize=9)
                plt.show()
            return

    def solveMatrix(self, cb):
        pivot_el_row = self.obj_func + self.real_vars + \
            self.pivot_el_col + self.constraint
        flag = 0
        sigma = -cb*self.Binv*self.A
        for i in range(len(self.pol)+1):
            for i in range(np.shape(self.sigma)[1]):
                if sigma[0, i] < 0:
                    flag = 1
                    break
            if not flag:
                break
            min_col_ind = np.where(sigma == np.amin(sigma))[1][0]
            pivot_col_val = self.A[:, min_col_ind]
            for i in range(len(self.A[:, min_col_ind])):
                if i == 0:
                    min_row_ind = i
                    if pivot_col_val[i] > 0:
                        b_a = self.xb[i, 0]/pivot_col_val[i]
                    else:
                        b_a = 100000
                else:
                    if pivot_col_val[i] > 0:
                        if self.xb[i, 0]/pivot_col_val[i] < b_a:
                            b_a = self.xb[i, 0]/pivot_col_val[i]
                            min_row_ind = i

            self.pivot_el_col[min_row_ind] = pivot_el_row[min_col_ind+1]
            cb[:, min_row_ind] = self.cd[:, min_col_ind]
            self.B[:, min_row_ind] = self.D[:, min_col_ind]
            self.Binv = np.linalg.inv(self.B)
            self.xb = self.Binv*self.b.T
            z = cb*self.Binv*self.b.T
            self.A = self.Binv*self.D

            for i in range(np.shape(self.A)[1]):
                suma = 0
                for j in range(np.shape(self.A)[0]):
                    if 'a' in self.pivot_el_col[j]:
                        suma += self.cd[0, i] - self.A[j, i]*1000*cb[0, j]
                    else:
                        suma += self.cd[0, i] - self.A[j, i]*cb[0, j]
                sigma[0, i] = suma

            res_dict = {}
            for i in range(len(self.pivot_el_col)-1):
                res_dict[self.pivot_el_col[i]] = np.asarray(
                    self.Binv*self.b.T)[i, :][0]

        return z
