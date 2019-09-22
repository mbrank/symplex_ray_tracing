class ConvexPolyhedron(object):
    """ConvexPolyhedron class

    """

    def __init__(self, pts):
        super(ConvexPolyhedron, self).__init__()
        self.pts = pts

    def __showPointsInfo__(self):
        for i in range(len(self.pts)):
            print("Index:", i, "point:",
                  self.pts[i][0],
                  self.pts[i][1])

    def addPoint(self, point):
        self.pts.append(point)

    def removePoint(self, index):
        self.pts.pop(index)

    def getListofPoints(self):
        return self.pts
