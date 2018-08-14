import math

def cleanPointList(list_of_points):
    """
    Deletes duplicate points in list_of_points
    """
    return list(set(list_of_points))


def addPointToList(vector, element):
    """
    Returns vector with the given element append to the right
    """
    vector.append(element)
    return vector


def removePointFromList(vector, element):
    """
    Returns a copy of vector without the given element
    """
    vector.pop(vector.index(element))
    return vector


def findMinYPoint(list_of_points):
    """
    Returns that point of *list_of_points* having minimal y-coordinate

    :param list_of_points: list of tuples
    :return: tuple (x, y)
    """
    min_y_pt = list_of_points[0]
    for point in list_of_points[1:]:
        if point[1] < min_y_pt[1] or (point[1] == min_y_pt[1] and point[0] < min_y_pt[0]):
            min_y_pt = point
    return min_y_pt

def calculateEuclideanDistance(point1, point2):
    """
    Returns the euclidian distance of the 2 given points.

    :param point1: tuple (x, y)
    :param point2: tuple (x, y)
    :return: float
    """
    return math.sqrt(math.pow(point1[0] - point2[0], 2) + math.pow(point1[1] - point2[1], 2))


def findNearestPoints(list_of_points, point, k):
    """
    gibt eine Liste mit den Indizes der k nächsten Nachbarn aus list_of_points zum angegebenen Punkt zurück.
    Das Maß für die Nähe ist die euklidische Distanz. Intern wird k auf das Minimum zwischen dem gegebenen Wert
    für k und der Anzahl der Punkte in list_of_points gesetzt

    :param list_of_points: list of tuples
    :param point: tuple (x, y)
    :param k: integer
    :return: list of k tuples
    """
    # build a list of tuples of distances between point *point* and every point in *list_of_points*, and
    # their respective index of list *list_of_distances*
    list_of_distances = []
    for index in range(len(list_of_points)):
        list_of_distances.append(( calculateEuclideanDistance(list_of_points[index], point), index))

    # sort distances in ascending order
    list_of_distances.sort()

    # get the k nearest neighbors of point
    nearest_list = []
    for index in range(min(k, len(list_of_points))):
        nearest_list.append((list_of_points[list_of_distances[index][1]]))
    return nearest_list


def sortByAngle(list_of_points, last_point, last_angle):
    """
    gibt die Punkte in list_of_points in absteigender Reihenfolge des Winkels zum letzten Segment der Hülle zurück,
    gemessen im Uhrzeigersinn. Es wird also immer der rechteste der benachbarten  Punkte ausgewählt. Der erste
    Punkt dieser Liste wird der nächste Punkt der Hülle.
    """
    def getkey(item):
        return calculateAngleDifference(last_angle, angle(last_point, item))

    vertex_list = sorted(list_of_points, key=getkey, reverse=True)
    return vertex_list


def intersect(line1, line2):
    """
    Returns True if the two given line segments intersect each other, and False otherwise.

    :param line1: 2-tuple of tuple (x, y)
    :param line2: 2-tuple of tuple (x, y)
    :return: boolean
    """
    a1 = line1[1][1] - line1[0][1]
    b1 = line1[0][0] - line1[1][0]
    c1 = a1 * line1[0][0] + b1 * line1[0][1]
    a2 = line2[1][1] - line2[0][1]
    b2 = line2[0][0] - line2[1][0]
    c2 = a2 * line2[0][0] + b2 * line2[0][1]
    tmp = (a1 * b2 - a2 * b1)
    if tmp == 0:
        return False
    sx = (c1 * b2 - c2 * b1) / tmp
    if (sx > line1[0][0] and sx > line1[1][0]) or (sx > line2[0][0] and sx > line2[1][0]) or\
            (sx < line1[0][0] and sx < line1[1][0]) or (sx < line2[0][0] and sx < line2[1][0]):
        return False
    sy = (a1 * c2 - a2 * c1) / tmp
    if (sy > line1[0][1] and sy > line1[1][1]) or (sy > line2[0][1] and sy > line2[1][1]) or\
            (sy < line1[0][1] and sy < line1[1][1]) or (sy < line2[0][1] and sy < line2[1][1]):
        return False
    return True

def angle(from_point, to_point):
    """
    Returns the angle of the directed line segment, going from *from_point* to *to_point*, in radians. The angle is
    positive for segments with upward direction (north), otherwise negative (south). Values ranges from 0 at the
    right (east) to pi at the left side (west).

    :param from_point: tuple (x, y)
    :param to_point: tuple (x, y)
    :return: float
    """
    return math.atan2(to_point[1] - from_point[1], to_point[0] - from_point[0])


def pointInPolygon(point, list_of_points):
    """
    Return True if given point *point* is laying in the polygon described by the vertices *list_of_points*,
    otherwise False

    Based on the "Ray Casting Method" described by Joel Lawhead in this blog article:
    http://geospatialpython.com/2011/01/point-in-polygon.html

    """
    x = point[0]
    y = point[1]
    poly = [(pt[0], pt[1]) for pt in list_of_points]
    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

def calculateAngleDifference(angle1, angle2):
    """
    Calculates the difference between the given angles in clockwise direction as radians.

    :param angle1: float
    :param angle2: float
    :return: float; between 0 and 2*Pi
    """
    if (angle1 > 0 and angle2 >= 0) and angle1 > angle2:
        return abs(angle1 - angle2)
    elif (angle1 >= 0 and angle2 > 0) and angle1 < angle2:
        return 2 * math.pi + angle1 - angle2
    elif (angle1 < 0 and angle2 <= 0) and angle1 < angle2:
        return 2 * math.pi + angle1 + abs(angle2)
    elif (angle1 <= 0 and angle2 < 0) and angle1 > angle2:
        return abs(angle1 - angle2)
    elif angle1 <= 0 < angle2:
        return 2 * math.pi + angle1 - angle2
    elif angle1 >= 0 >= angle2:
        return angle1 + abs(angle2)
    else:
        return 0


def calcualateAlphaShape(points_list, k):
    """
    Calculates a valid concave hull polygon containing all given points. The algorithm searches for that
    point in the neighborhood of k nearest neighbors which maximizes the rotation angle in clockwise direction
    without intersecting any previous line segments.

    This is an implementation of the algorithm described by Adriano Moreira and Maribel Yasmina Santos:
    CONCAVE HULL: A K-NEAREST NEIGHBOURS APPROACH FOR THE COMPUTATION OF THE REGION OCCUPIED BY A SET OF POINTS.
    GRAPP 2007 - International Conference on Computer Graphics Theory and Applications; pp 61-68.

    :param points_list: list of tuples (x, y)
    :param k: integer
    :return: list of tuples (x, y)
    """
    # return an empty list if not enough points are given
    if k > len(points_list):
        return None

    # the number of nearest neighbors k must be greater than or equal to 3
    # kk = max(k, 3)
    kk = max(k, 2)

    # delete duplicate points
    point_set = cleanPointList(points_list)

    # if point_set has less then 3 points no polygon can be created and an empty list will be returned
    if len(point_set) < 3:
        return None

    # if point_set has 3 points then these are already vertices of the hull. Append the first point to
    # close the hull polygon
    if len(point_set) == 3:
        return addPointToList(point_set, point_set[0])

    # make sure that k neighbours can be found
    kk = min(kk, len(point_set))

    # start with the point having the smallest y-coordinate (most southern point)
    first_point = findMinYPoint(point_set)

    # add this points as the first vertex of the hull
    hull = [first_point]

    # make the first vertex of the hull to the current point
    current_point = first_point

    # remove the point from the point_set, to prevent it being among the nearest points
    point_set = removePointFromList(point_set, first_point)
    previous_angle = math.pi

    # step counts the number of segments
    step = 2

    # as long as point_set is not empty or search is returning to the starting point
    while (current_point != first_point) or (step == 2) and (len(point_set) > 0):

        # after 3 iterations add the first point to point_set again, otherwise a hull cannot be closed
        if step == 5:
            point_set = addPointToList(point_set, first_point)

        # search the k nearest neighbors of the current point
        k_findNearestPoints = findNearestPoints(point_set, current_point, kk)

        # sort the candidates (neighbors) in descending order of right-hand turn. This way the algorithm progresses
        # in clockwise direction through as many points as possible
        c_points = sortByAngle(k_findNearestPoints, current_point, previous_angle)

        its = True
        i = -1

        # search for the nearest point to which the connecting line does not intersect any existing segment
        while its is True and (i < len(c_points)-1):
            i += 1
            if c_points[i] == first_point:
                last_point = 1
            else:
                last_point = 0
            j = 2
            its = False

            while its is False and (j < len(hull) - last_point):
                its = intersect((hull[step-2], c_points[i]), (hull[step-2-j], hull[step-1-j]))
                j += 1

        # there is no candidate to which the connecting line does not intersect any existing segment, so the
        # for the next candidate fails. The algorithm starts again with an increased number of neighbors
        if its is True:
            return calcualateAlphaShape(points_list, kk + 1)

        # the first point which complies with the requirements is added to the hull and gets the current point
        current_point = c_points[i]
        hull = addPointToList(hull, current_point)

        # calculate the angle between the last vertex and his precursor, that is the last segment of the hull
        # in reversed direction
        previous_angle = angle(hull[step - 1], hull[step - 2])

        # remove current_point from point_set
        point_set = removePointFromList(point_set, current_point)

        # increment counter
        step += 1

    all_inside = True
    i = len(point_set)-1

    # check if all points are within the created polygon
    while (all_inside is True) and (i >= 0):
        all_inside = pointInPolygon(point_set[i], hull)
        i -= 1

    # if at least one point is out of the computed polygon, try again with a higher number of neighbors
    if all_inside is False:
        return calcualateAlphaShape(points_list, kk + 1)

    # a valid hull has been constructed
    return hull