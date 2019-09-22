from ComputeCrossSection import ComputeCrossSection
from Ray import Ray
from ConvexPolyhedron import ConvexPolyhedron
import numpy as np


def main():
    ray = Ray(np.matrix([-1, -1]), np.matrix([1, 1]))
    pol = ConvexPolyhedron(np.matrix([[-10, 20, 50, 30, 10],
                                      [20, 40, 30, 10, -10]]))
    cs = ComputeCrossSection(ray, pol)
    cs.checkIntersection()


if __name__ == '__main__':
    main()
