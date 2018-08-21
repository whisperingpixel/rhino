import gdal
import numpy
import shapely
import shapely.geometry as shapely_geometry
import shapely.ops as operations
import math
import copy

from rhino_helper import progressBar
from rhino_helper import checkCoordinates

## TODO:
# - generate global id for atoms
# - createObjectView function is misleading
# - Handle different bands
# - Handle time

class Datacube:

    dataset = None

    config = {}

    object_model = {
        "derived":{
            "compactness": None,
            "boundingbox": None,
            "geometry": None
        },
        "modelled": {
            "class": None
        }
    }

    link_model = {
        "derived": ["link_type"],
        "modelled": []
    }

    neighbourhood = None

    def __init__(self):

        self.__view = {
            "coverage": [],
            "objects": []
        }

        self.__levels = []

        #
        # TODO: Make proper config
        #
        Datacube.config["showProgress"] = False

    def load(self, connection, boundingbox):
        """
        Loads the rhino datacube.

        :param boundingbox: [x_min, x_max, y_min, y_max]
        """
        ## TODO: Load datacube from here. Be advised that this is a workaraound without a datacube

        Datacube.dataset = gdal.Open(connection)
        
        ## TODO: account for different bands

        level = self.createNewLevel("initial_level")
        #
        # Note there are some special treatments since this is the loading procedure (hen-egg problem). Do not 
        # do this somewhere else. Basically it avoids that the atoms have to be created multiple times. The numpy array 
        # should not be further exposed to the user as the theory says it is actually a field of atoms.
        #
        band = Datacube.dataset.GetRasterBand(1)
        initial_coverage = Coverage(Datacube.dataset.RasterXSize, Datacube.dataset.RasterYSize, Datacube.dataset.GetGeoTransform())
        initial_coverage.create(array = numpy.array(band.ReadAsArray()))
        self.__view["coverage"][level] = initial_coverage
        self.createObjectView(level, level)

        #
        # Use the atoms created by the object view for now
        #
        atoms = []
        for o in self.__view["objects"][0]:
            atoms.append(o.getAtoms()[0]) #there is only one atom per object now
        self.__view["coverage"][0].create(atoms = atoms)

        #
        # Set default neighbourhood concept
        #
        Datacube.neighbourhood = Neighbourhood("4-connected")


    def createNewLevel(self, name ="no name available", description = "no description available"):
        """
        Creates a new level, which is pre-condition for creating any object or coverage view.

        :param name: string; name of the level (optional)
        :param description: string; decription of the level (optional)
        """        
        new_level = len(self.__levels)
        self.__levels.append({
            "level": new_level,
            "name": name,
            "description": description 
        })

        self.__view["coverage"].append(Coverage)
        self.__view["objects"].append([])

        return new_level


    def getLevel(self, depth):
        """
        Returns the level by the depth

        :param depth: integer
        :return: dict
        """
        for level in self.__levels:
            if level["level"] == depth:
                return level
        return None


    def createObjectView(self, from_level, to_level):
        """
        Creates the object view in the rhino datacube.

        :param layer: int
        """

        #
        # Init an empty object view
        #
        object_list = []

        # 
        # Get metadata (esp. getransform) from the coverage
        #
        (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = self.__view["coverage"][from_level].getMetadata()

        #
        # Iterate through all pixels and create an atom from it. Then, create an object for each atom.
        # Initial objects only consist of one single atom
        #
        # Stolen from here: https://stackoverflow.com/questions/6967463/iterating-over-a-numpy-array/6967491#6967491

        if Datacube.config["showProgress"]:
            metadata = self.__view["coverage"][from_level].getSize()
            size = metadata[0] * metadata[1]
            counter = 0            
            progressBar(counter, size, prefix = 'Generate views:', suffix = 'Complete', length = 50)
        
        for (x_index, y_index), value in numpy.ndenumerate(self.__view["coverage"][from_level].getStupidArray()):
            
            #
            # Extract coordinates, properties and values from atoms
            # Stolen from here: https://gis.stackexchange.com/questions/42790/gdal-and-python-how-to-get-coordinates-for-all-cells-having-a-specific-value/42846#42846

            #
            # TODO: There is an issue with x/y lat/lon here.
            #
            lon = y_index * x_size + upper_left_x + (x_size / 2) #add half the cell size
            lat = x_index * y_size + upper_left_y + (y_size / 2) #to centre the point
            # TODO: add time
            property = "class"
            value = int(value) ## TODO: remove int, this is due to the demo data

            #
            # Create an object and fill it with the atom
            #
            o = Object()
            o.create([Atom({"lat": lat, "lon": lon},{"property": property, "value": value},{"x":x_index,"y":y_index})])
            
            object_list.append(o)

            if Datacube.config["showProgress"]: 
                counter += 1                          
                progressBar(counter, size, prefix = 'Generate views:', suffix = 'Complete', length = 50)
    
        self.__view["objects"][to_level] = object_list


    def createCoverageView(self, from_level, to_level):
        """
        Creates the coverage view in the rhino datacube.
        This can also be seen as masking, e.g., after selection of objects.

        :param objects: list of objects
        """
        old_stupid_array = self.__view["coverage"][0].getStupidArray()
        x_size = len(old_stupid_array)
        y_size = len(old_stupid_array[0])

        new_stupid_array = numpy.array(old_stupid_array, copy=True)
        new_stupid_array = new_stupid_array * 0
        
        atoms = []
        for o in self.__view["objects"][from_level]:
            for a in o.getAtoms():
                x,y = a.getIndex()
                value = a.getObservationValue()
                new_stupid_array[x][y] = value
                atoms.append(a)

        coverage = Coverage(x_size, y_size, self.__view["coverage"][0].getMetadata())
        coverage.create(atoms=atoms, array=new_stupid_array)

        self.__view["coverage"][to_level] = coverage


    def getDatasetMetadata(self):
        """
        Returns the metadata of the selected dataset.

        :return: dict
        """
        gt = Datacube.dataset.GetGeoTransform()
        metadata = {
            "coverage_x_size": Datacube.dataset.RasterXSize,
            "coverage_y_size": Datacube.dataset.RasterYSize,
            "coverage_x_grain": gt[1],
            "coverage_y_grain": gt[5]
        }

        return metadata

    ## TODO: This is currently not possible as switching to another neighbourhood
    #        would invalidate existing object links.
    #def setNeighbourhoodConcept(self, concept):
    #    neighbourhood = Neighbourhood(concept)
    

    #
    # These are user-friendly wrapper-functions
    #
    def selectObjectsByCondition(self, from_level, condition, create_level = False):
        """
        Select the objects in the object view of the rhino datacube, which match
        a certain condition.

        :param from_level: integer, the level which should be used to select the objects
        :param condition: ?
        :param create_level: boolean, whether a new level should be created or not
        :return: list of objects
        """

        def eq(value, candidate):
            return value == candidate

        def neq(value, candidate):
            return value != candidate

        def lt(value, candidate):
            return value < candidate

        def let(value, candidate):
            return value <= candidate

        def gt(value, candidate):
            return value > candidate

        def get(value, candidate):
            return value >= candidate

        objects = []
        candidates = self.__view["objects"][from_level]

        if Datacube.config["showProgress"]:
            counter = 0
            size = len(candidates)
            progressBar(counter, size, prefix = 'Selecting objects with condition:', suffix = 'Complete', length = 50)

        for obj in candidates:
            if locals()[condition["operator"]](obj.getAttributeValue(condition["key"]), condition["value"]): ## TODO: from condition!
                objects.append(obj)
            
            if Datacube.config["showProgress"]:
                counter += 1
                progressBar(counter, size, prefix = 'Selecting objects with condition:', suffix = 'Complete', length = 50)


        if create_level == True:

            level = self.createNewLevel()
            self.__view["objects"][level] = objects
            self.createCoverageView(level, level)
            coverage = self.__view["coverage"][level]

        else:

            level = None
            coverage = None
    
        return objects, coverage, level


    def aggregate(self, from_level, condition):
        """
        Aggregates ? based on a certain condition.

        :param condition: ?
        :return: ?
        """

        ## TODO: use create new level swthc in select method
        #
        # (optional): Create new level. In this case, we use the interal method in the next step.
        #

        #
        # STEP 1: Select all objects, which match the condition. If no objects
        # match the condition, return here.
        #
        # This function generates a new level and views
        #
        objects, dummy_coverage, level = self.selectObjectsByCondition(from_level, condition, create_level=True)
        if len(objects) == 0:
            return

        #
        # STEP 2: Aggregate and update views
        #

        #
        # Aggregate/cluster using the Link class, it updates the levels accordingly
        #
        Link("aggregation", self.__view, level)

        #
        # (Optional) return information about the new level
        #
        return self.getLevel(level)


    def executeQuery(self, command):
        """
        Executes a SQL-like query command against the object view of the rhino data cube.

        :param command: string
        :return: ?
        """
        pass ##TODO: Implement


class Atom:

    def __init__(self, coordinates, observation, index):

        ## TODO: Check value integrity
        checkCoordinates(coordinates["lat"], coordinates["lon"])

        self.__tuple = {
            "id": id(self),
            "lat": coordinates["lat"],
            "lon": coordinates["lon"],
            "geometry": shapely_geometry.Point(coordinates["lat"],coordinates["lon"]),
            "property": observation["property"],
            "value": observation["value"],
            "x_index": index["x"],
            "y_index": index["y"]
        }


    def getID(self):
        """
        Returns the id of the atom.

        :return: int
        """        
        return self.__tuple["id"]


    def getObservationProperty(self):
        """
        Returns the observation property of the atom.

        :return: string
        """
        return self.__tuple["property"]


    def getObservationValue(self):
        """
        Returns the observation value of the atom.

        :return: ?
        """        
        return self.__tuple["value"]


    def getObservation(self):
        """
        Returns the observation of the atom as dictionary:
        {"property": property, "value": value}

        :return: dict
        """               
        return {"property": self.__tuple["property"], "value": self.__tuple["value"]}


    def getLatitude(self):
        """
        Returns the latitude of the atom.

        :return: float
        """        
        return self.__tuple["lat"]


    def getLongitude(self):
        """
        Returns the longitude of the atom.

        :return: float
        """
        return self.__tuple["lon"]


    def getCoordinates(self):
        """
        Returns the coordinates of the atom as tuple

        :return: tuple(lat, lon)
        """        
        return (self.__tuple["lat"], self.__tuple["lon"])


    def getGeometry(self):
        """
        Returns the geometry of the atom as shapely point

        :return: shapely point
        """
        return self.__tuple["geometry"]

    def getIndex(self):
        """
        Returns the internal index (position within the coverage) for that atom

        :return: shapely point
        """
        return (self.__tuple["x_index"], self.__tuple["y_index"])

    def print(self):
        print(self.__tuple)

class Coverage:

    def __init__(self, x_size, y_size, geotransform):
        self.id = id(self)
        self.coverage = [[]]
        self.__numpyarray = None
        self.geotransform = geotransform

        self.x_size = x_size
        self.y_size = y_size

        for x in range(x_size):
            self.coverage.append([])
            for y in range(y_size):
                self.coverage[x].append(None)
    
    def create(self, atoms = None, array = None):
        """
        Creates the coverage (create as in "creation" by god), don't confuse with
        initialisation. It takes a list of atoms and puts it into a coverage.

        One of the following parameters has to be set.
        :param atoms: list of atoms
        :param array: numpy array
        """
        if atoms is None and array is None:
            raise Exception("Coverage creation failed")

        if atoms is not None:
            for atom in atoms:
                x,y = atom.getIndex()
                self.coverage[x][y] = atom

        if array is not None:
            self.__numpyarray = array

    def getCoverage(self, boundingbox=None):
        """
        Returns the atom coverage.

        :param boundingbox: ?
        :return: atom coverage
        """        
        if boundingbox is None:
            return self.coverage
        else:
            return self.coverage
    
    def getStupidArray(self):
        """
        Returns the stupid (numpy) array. This should only be used for
        performance reasons and hidden from the user.

        :return: numpy array
        """         
        return self.__numpyarray

    def getSize(self):
        """
        Returns the size in x and y direction of the coverage.

        :return: tuple (x_size, y_size)
        """                 
        return (self.x_size, self.y_size)


    def getMetadata(self):
        """
        Returns the metadata (geotransform) of the coverage.

        :return: tuple (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size)
        """                 
        return self.geotransform

class Object:

    def __init__(self):
        self.id = id(self)
        self.atoms = []
        self.links = []
        self.__attributes = copy.deepcopy(Datacube.object_model)

    def create(self, atoms):
        """
        Creates the object (create as in "creation" by god), don't confuse with
        initialisation. It takes a list of atoms and derives geospatial properties of the object.

        :param atoms: list of atoms
        """
        self.atoms = atoms
#        self.__calculateAlphaShape(3)
        self.__calculateGeoAttributes()

    def grow(self, atoms):
        """
        Grows the object with the given atomes

        :param atoms: list of atoms
        """
        self.atoms.extend(atoms)
#        self.__calculateAlphaShape(3)
        self.__calculateGeoAttributes()
        self.__attributes["derived"]["geometry"] = None

    def shrink(self, atoms):
        """
        Shrinks the object by the given atomes

        :param atoms: list of atoms
        """
        self.__attributes["derived"]["geometry"] = None

    def __calculateAlphaShape(self, neighbors):
        
        from rhino_geo import calcualateAlphaShape
        polygon = calcualateAlphaShape(self.getAtomCoordinates(), neighbors)
        if polygon:
            self.__attributes["derived"]["outline"] = polygon
        ## TODO: Handle case where not enough atoms are available for that polygon

    def __calculateGeoAttributes(self):

        self.__attributes["derived"]["number_of_atoms"] = len(self.atoms)

    def getAtoms(self):
        """
        Returns the atoms, which are forming the object.

        :return: list of atoms
        """
        return self.atoms
    
    def getNumberOfAtoms(self):
        """
        Returns the number of atoms, which are forming the object.

        :return: integer
        """
        return len(self.atoms)

    def getAtomCoordinates(self):
        """
        Returns the coordinates of the atoms, which are forming the object, as list.

        :return: list of tuples (lat, lon)
        """
        coordinates = []
        for a in self.atoms:
            coordinates.append(a.getCoordinates())
        return coordinates

    def getBoundingBox(self):
        """
        Returns bounding box of the object as shapely polygon

        :return: shapely polygon
        """
        if self.__attributes["derived"]["boundingbox"] == None:
            points = []
            for a in self.atoms:
                points.append(a.getGeometry())
            multipoints = shapely_geometry.MultiPoint(list(points))
            boundingbox = multipoints.envelope
            self.__attributes["derived"]["boundingbox"] = boundingbox
        else:
            boundingbox = self.__attributes["derived"]["boundingbox"]

        return boundingbox
    
    def getGeometry(self):
        """
        Returns geometry of the object as shapely polygon. Do not use this, it is very slow

        :return: shapely polygon
        """
        if self.__attributes["derived"]["geometry"] == None:
            points = []
            for a in self.atoms:
                points.append(a.getGeometry())
            multipoints = shapely_geometry.MultiPoint(list(points))

            #
            # TODO: Search for better algorithm
            #
            buf = Datacube.neighbourhood.getMaxDistance()
            geometry = operations.cascaded_union(multipoints.buffer(buf))
            geometry = geometry.buffer(- (buf / 2))
            self.__attributes["derived"]["geometry"] = geometry
        else:
            geometry = self.__attributes["derived"]["geometry"]

        return geometry

    def linksWithObjects(self, link_type, candidate):
        """
        Creates a link to an object candidate.

        :param link_type: string, type of the link: One of ["contract", "constraint","generalisation"]
        """
        def createLink(link_type, candidate):
            self.links.append({
                "link_type": link_type,
                "targets": candidate
            })

        if len(self.links) == 0:
            createLink(link_type, candidate)
        else:
            for link in self.links:
                if link["link_type"] == link_type:
                    link["targets"].append(candidate)
                else:
                    createLink(link_type, candidate)

    def getID(self):
        """
        Returns the id of the object.

        :return: int
        """
        return self.id

    def getAttributeValue(self, attribute):
        if attribute in self.__attributes["derived"]:
            if attribute == "compactness":

                geom = self.getGeometry()
                return (4 * math.pi * geom.area) / (geom.length * geom.length)

        else:
            if attribute in self.__attributes["modelled"]:
                atomProperties = []
                for a in self.atoms:
                    atomProperties.append(a.getObservationValue()) ## TODO: return based on the attribute arg
                return atomProperties[0] ## TODO: this returns the first one, make it user-defined                
            else:
                return None
        

    def setAttributeValue(self, attribute, value):
        """
        Sets the value of an attribute of this object. The attribute has to be in the
        object model for modelled attributes.

        :param attribute: string
        :param value: ?
        """        
        if attribute in self.__attributes["modelled"]:
            self.__attributes["modelled"][attribute] = value
        ## TODO: rais exception


    def print(self):
        """
        Prints the object.
        """                
        print("{ id: " + str(self.id) + ", atoms: ") 
        for a in self.atoms:
            a.print()
        print(",links: ")
        for l in self.links:
            print(l)
        print("}")


class Link:

    def __init__(self, link_type, target, level, mutual = True):

        self.id = id(self)
        self.attributes = {"derived":{}}
        self.attributes["derived"]["link_type"] = link_type

        self.__new_objects = []

        if link_type == "aggregation":

            self.__aggregation(target["coverage"][level])
            target["objects"][level] = self.__new_objects
        
        else:
            pass

#   Do stuff here
#            target[0].linksWithObjects(link_type, target[1])
#            if mutual == True:
#                target[1].linksWithObjects(link_type, target[0])


    def __aggregation(self, coverage):
        """
        Aggregates the atoms to objects. Background is considered 0!

        :param coverage: Coverage
        """

        stupid_array = coverage.getStupidArray()
        atom_coverage = coverage.getCoverage(None)

        x_size, y_size = coverage.getSize()

        #
        # Get regions using the stupid array
        #
        from skimage.measure import label as region_group

        grouped = region_group(stupid_array)

        #
        # Create a nested list where the atoms will be stored
        #
        atom_lists = []
        for i in range(numpy.max(grouped)):
            atom_lists.append([])

        #
        # Aggregate the atoms to the atom_list
        #
        if Datacube.config["showProgress"]:
            counter = 0
            progressBar(counter, x_size*y_size, prefix = 'Object linking (aggregate)', suffix = 'Complete', length = 50)

        for x in range(x_size):
             for y in range(y_size):
                atom = atom_coverage[x][y]
                object_number = grouped[x][y]
                if object_number != 0:
                    atom_lists[object_number-1].append(atom)

                if Datacube.config["showProgress"]:
                    counter += 1
                    progressBar(counter, x_size*y_size, prefix = 'Object linking (aggregate)', suffix = 'Complete', length = 50)
            
        #
        # Make objects
        #
        for atoms in atom_lists:
            o = Object()
            o.create(atoms)
            self.__new_objects.append(o)


    def getID(self):
        """
        Returns the id of the link (relationship/aggregation)

        :return: int
        """                
        return self.id

class Neighbourhood:

    def __init__(self, concept):

        if concept in ["raster_4", "4-connected", "raster_8", "moore", "8-connected"]:
            self.concept = concept
        else:
            raise Exception("Neighborhood concept " + str(concept) + " not implemented")
        
        self.tolerance = 0.001
    
    def getConcept(self):
        """
        Returns selected concept of neighbourhood

        :return: string
        """                        
        return self.concept

    def getMaxDistance(self):

        ## TODO: Account for different grains in X and Y direction
        dist = max(Datacube.getDatasetMetadata(Datacube)["coverage_x_grain"], Datacube.getDatasetMetadata(Datacube)["coverage_y_grain"])

        #
        # add tolerance
        #
        dist += self.tolerance

        return dist

    def isNeighbourObject(self, target, candidate):
        """
        Compares two objects and evaluates whether they are neighbours or not, based on the
        selected concept of neighbourhood.

        :param target: Target object
        :param candidate: 
        :return: boolean
        """                
        min_dist = float("inf")

        #
        # Get allowed distance
        #
        allowed_dist = self.getMaxDistance()

        #
        # Extract the atoms
        #
        obj_atoms = target.getAtoms()
        cand_atoms = candidate.getAtoms()

        #
        # Check the bounding boxes at first. If their distance is larger than the allowed distance
        # it is impossible that the object's distance is smaller than the allowed distance
        #
        if target.getBoundingBox().distance(candidate.getBoundingBox()) < allowed_dist:
            for a1 in obj_atoms:
                for a2 in cand_atoms:
                    dist = a1.getGeometry().distance(a2.getGeometry())
                    if dist < min_dist:
                        min_dist = dist
        else:
            return False

        #
        # Check whether the detected distance is smaller than the allowed distance
        #
        if min_dist < allowed_dist: ##TODO: this is for testing
            return True
        else:
            return False

    def isNeighbourInCoverage(self, target, candidate):
        pass