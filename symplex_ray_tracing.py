import numpy as np
import fractions
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
# np.set_printoptions(formatter={'all': lambda x: str(fractions.Fraction(x).limit_denominator())})



###############################################################################
#                               Define polygons                               #
###############################################################################

pol1 = np.matrix([[-10, 20, 50, 30, 10],
                  [20, 40, 30, 10, -10]])

pol2 = np.matrix([[40, 50, 10, -10],
                 [-10, 10, 40, 50]])

pol3 = np.matrix([[70, 50, 0],
                 [-10, 20, 40]])

pol4 = np.matrix([[-5, -5, 25],
                 [9.99999, 30, 0]])

pol_list = [pol1, pol2, pol3, pol4]

###############################################################################
#                                 Define rays                                 #
###############################################################################

Ray1 = np.matrix([1, 1])
Ray01 = np.matrix([-1, -1])

Ray2 = np.matrix([60, 10])
Ray02 = np.matrix([-1, -1])

Ray3 = np.matrix([50, 10])
Ray03 = np.matrix([-1, -1])

Ray4 = np.matrix([30, 10])
Ray04 = np.matrix([-1, -1])

Ray5 = np.matrix([2, 1])
Ray05 = np.matrix([-1, -1])

Ray6 = np.matrix([7, 1])
Ray06 = np.matrix([-1, -1])

Ray7 = np.matrix([8, 1])
Ray07 = np.matrix([-1, -1])

Ray8 = np.matrix([2, -1])
Ray08 = np.matrix([-1, -1])

Ray9 = np.matrix([80, 1])
Ray09 = np.matrix([-1, -1])

Ray10 = np.matrix([9, 1])
Ray010 = np.matrix([-1, -1])

Ray11 = np.matrix([10, 1])
Ray011 = np.matrix([-1, -1])

Ray12 = np.matrix([11, 1])
Ray012 = np.matrix([-1, -1])

Ray13 = np.matrix([12, 1])
Ray013 = np.matrix([-1, -1])

Ray14 = np.matrix([14, 1])
Ray014 = np.matrix([-1, -1])

Ray15 = np.matrix([1, 10])
Ray015 = np.matrix([-1, -1])

Ray16 = np.matrix([1, 15])
Ray016 = np.matrix([-1, -1])

Ray17 = np.matrix([1, 11])
Ray017 = np.matrix([-1, -1])

Ray18 = np.matrix([4, 11])
Ray018 = np.matrix([-1, -1])

Ray19 = np.matrix([-10, 10])
Ray019 = np.matrix([-1, -1])

Ray20 = np.matrix([-10, -10])
Ray020 = np.matrix([-1, -1])

Ray_list = [Ray1, Ray2, Ray3, Ray4, Ray5, Ray6, Ray7, Ray8, Ray9,
            Ray10, Ray11, Ray12, Ray13, Ray14, Ray15, Ray16, Ray17, Ray18,
            Ray19, Ray20]
Ray0_list = [Ray01, Ray02, Ray03, Ray04, Ray05, Ray06, Ray07, Ray08, Ray09,
             Ray010, Ray011, Ray012, Ray013, Ray014, Ray015, Ray016, Ray17,
             Ray018, Ray019, Ray020]

fig, ax = plt.subplots()
patches = []
num_polygons = 5
num_sides = 5

###############################################################################
#                  Iterate through polygons and find intersection             #
###############################################################################

for l in range(len(pol_list)):
    pol = pol_list[l]

    polygon = Polygon(pol.T, True)
    patches.append(polygon)
    p = PatchCollection(patches, alpha=0.4)
    ax.add_collection(p)
    plt.xlim(-20, 80)
    plt.ylim(-20, 50)

    for pt in range(np.shape(pol)[1]):
        plt.scatter(pol[0, pt], pol[1, pt])
        plt.text(pol[0, pt]+0.1, pol[1, pt]+0.1, "P"+str(pt-1), fontsize=9)


    # iteration through rays
    for k in range(len(Ray_list)):
        Ray = Ray_list[k]
        Ray0 = Ray0_list[k]

        plt.scatter(Ray[0, 0]+Ray0[0, 0], Ray[0, 1]+Ray0[0, 1])
        plt.text(Ray[0, 0]+Ray0[0, 0]+0.1, Ray[0, 1]+Ray0[0, 1]+0.1, "R"+str(k+1), fontsize=9)
        plt.scatter(Ray0[0, 0], Ray0[0, 1])
        plt.plot([Ray0[0, 0], Ray0[0, 0]+Ray[0, 0]],
                 [Ray0[0, 1], Ray0[0, 1]+Ray[0, 1]])

        cd = np.block([np.matrix([0 for i in range(np.shape(pol)[1])]),
                       np.matrix([0])])

        cb = np.matrix([1, 1, 1])

        A = np.block([[pol, -Ray.T],
                      [np.matrix([1 for i in range(np.shape(pol)[1])]),
                       np.matrix([0])]])

        b = np.matrix([[Ray0[0, 0], Ray0[0, 1], 1]])

        B = np.matrix([[1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])
        Binv = np.linalg.inv(B)

        for i in range(np.shape(b)[1]-1):
            if b[0, i] < 0:
                b[0, i] *= -1
                A[i] *= -1

        D = A

        sigma = -cb*Binv*A

        xb = Binv*b.T
        z = cb*Binv*b.T

        obj_func = ['z']
        real_vars = ['x'+str(i) for i in range(np.shape(pol)[1]+1)]
        lam = real_vars[-1]
        art_vars = ['ax'+str(i) for i in range(np.size(b))]
        constraint = ['b']
        pivot_el_row = obj_func + real_vars + art_vars + constraint
        pivot_el_col = art_vars
        flag = 0
        for i in range(20):
            flag = 0
            for i in range(np.shape(sigma)[1]):
                if sigma[0, i] < 0:
                    flag = 1
                    break
            if not flag:
                break
            min_col = np.asarray(sigma)[0]
            min_col_ind = np.where(sigma == np.amin(sigma))[1][0]
            pivot_col_val = A[:, min_col_ind]
            for i in range(len(A[:, min_col_ind])):
                if i == 0:
                    min_row_ind = i
                    if pivot_col_val[i] > 0:
                        b_a = xb[i, 0]/pivot_col_val[i]
                    else:
                        b_a = 100000
                else:
                    if pivot_col_val[i] > 0:
                        if xb[i, 0]/pivot_col_val[i] < b_a:
                            b_a = xb[i, 0]/pivot_col_val[i]
                            min_row_ind = i

            pivot_el_col[min_row_ind] = pivot_el_row[min_col_ind+1]
            B[:, min_row_ind] = D[:, min_col_ind]
            cb[:, min_row_ind] = cd[:, min_col_ind]
            Binv = np.linalg.inv(B)
            xb = Binv*b.T
            z = cb*Binv*b.T
            A = Binv*D

            for i in range(np.shape(A)[1]):
                suma = 0
                for j in range(np.shape(A)[0]):
                    if 'a' in pivot_el_col[j]:
                        suma += cd[0, i] - A[j, i]*1000*cb[0, j]
                    else:
                        suma += cd[0, i] - A[j, i]*cb[0, j]
                sigma[0, i] = suma

            res_dict = {}
            for i in range(len(pivot_el_col)-1):
                res_dict[pivot_el_col[i]] = np.asarray(Binv*b.T)[i, :][0]

        if z[0, 0] != 0:
            print("Ray "+str(k+1) + " and Polygon "+str(l+1)+": "
                  "No intersection, Z NOT EQUAL TO 0")

        else:
            cd = np.block([np.matrix([0 for i in range(np.shape(pol)[1])]),
                           np.matrix([1])])

            cb = np.matrix([1, 1, 1])

            sigma = -cb*Binv*A
            flag = 0
            for i in range(4):
                flag = 0
                for i in range(np.shape(sigma)[1]):
                    if sigma[0, i] < 0:
                        flag = 1
                        break
                if not flag:
                    break
                min_col = np.asarray(sigma)[0]
                min_col_ind = np.where(sigma == np.amin(sigma))[1][0]
                pivot_col_val = A[:, min_col_ind]
                for i in range(len(A[:, min_col_ind])):
                    if i == 0:
                        min_row_ind = i
                        if pivot_col_val[i] > 0:
                            b_a = xb[i, 0]/pivot_col_val[i]
                        else:
                            b_a = 100000
                    else:
                        if pivot_col_val[i] > 0:
                            if xb[i, 0]/pivot_col_val[i] < b_a:
                                b_a = xb[i, 0]/pivot_col_val[i]
                                min_row_ind = i

                pivot_el_col[min_row_ind] = pivot_el_row[min_col_ind+1]
                B[:, min_row_ind] = D[:, min_col_ind]
                cb[:, min_row_ind] = cd[:, min_col_ind]
                Binv = np.linalg.inv(B)
                xb = Binv*b.T
                z = cb*Binv*b.T
                A = Binv*D

                for i in range(np.shape(A)[1]):
                    suma = 0
                    for j in range(np.shape(A)[0]):
                        if 'a' in pivot_el_col[j]:
                            suma += cd[0, i] - A[j, i]*1000*cb[0, j]
                        else:
                            suma += cd[0, i] - A[j, i]*cb[0, j]
                    sigma[0, i] = suma

            res_dict = {}
            for i in range(len(pivot_el_col)-1):
                res_dict[pivot_el_col[i]] = np.asarray(Binv*b.T)[i, :][0]

            if lam in res_dict:
                ex1 = Ray0[0, 0]+Ray[0, 0]*res_dict[lam]
                ey1 = Ray0[0, 1]+Ray[0, 1]*res_dict[lam]
                plt.scatter(ex1, ey1)
                plt.text(ex1+0.1, ey1+0.1, "e"+str(k+1), fontsize=9)
            else:
                pass


plt.grid()
plt.show()
