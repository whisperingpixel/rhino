import gdal
import numpy
import shapely
import shapely.geometry as shapely_geometry
import shapely.ops as shapely_operations
from shapely.wkt import loads as load_wkt
import math
import copy
import datacube
import xarray
import re

from . import rhino_tools

## TODO:
# - generate global id for atoms
# - Handle different bands
# - Handle time

class Datacube:

    dataset = None

    config = {}

    object_model = {
        "derived":{
            "compactness": None,
            "area": None,
            "boundingbox": None,
            "geometry": None
        },
        "modelled": {
            "ndvi": None
        }
    }

    link_model = {
        "derived": {
            "link_type": None
        },
        "modelled": []
    }
    
    dataset_metadata = {
        "coverage_x_size": None,
        "coverage_y_size": None,
        "coverage_x_grain": None,
        "coverage_y_grain": None
    }

    neighbourhood = None
    atom_property = None


    def __init__(self, show_progress = False, atom_property = None):

        self.__view = {
            "coverage": [],
            "objects": []
        }

        self.__levels = []

        #
        # TODO: Make proper config
        #
        Datacube.config["show_progress"] = show_progress
        Datacube.config["sfc_index"] = False

        if atom_property is not None:
            Datacube.atom_property = atom_property

        #
        # Indexing using Hilbert Curve
        #
        self.i_xmin = None
        self.i_xmax = None
        self.i_ymin = None
        self.i_ymax = None
        self.i_tmax = None
        self.i_tmin = None
        self.i_number_of_cells = 100
        self.sfc = None
        self.sfc_index = []

    def load(self, product, domain):
        """
        Loads the rhino datacube.

        :param domain: [x_min, x_max, y_min, y_max]
        """

        domain_shape = load_wkt(domain["geometry"])
        domain_bounds = domain_shape.bounds
        domain_crs = domain["crs"]

        dc = datacube.Datacube()
        area_of_interest = dc.load(
            product = product,
            x = (domain_bounds[0], domain_bounds[2]),
            y = (domain_bounds[1],domain_bounds[3]),
            crs = domain_crs)

        if Datacube.config["sfc_index"]:
            self.i_xmin = float(min(area_of_interest.coords["y"]))
            self.i_xmax = float(max(area_of_interest.coords["y"]))
            self.i_ymin = float(min(area_of_interest.coords["x"]))
            self.i_ymax = float(max(area_of_interest.coords["x"]))
            self.t_min = 0.0
            self.t_max = 1.0

            self.sfc = rhino_tools.HilbertCurve(self.i_number_of_cells, 3)
            for i in range(self.getIndexPosition((self.i_xmax, self.i_ymax, self.i_tmax))):
                self.sfc_index.append([])

        self.setDatasetMetadata(
            coverage_x_size = area_of_interest.sizes["x"],
            coverage_y_size = area_of_interest.sizes["y"],
            coverage_t_size = area_of_interest.sizes["time"],
            coverage_x_grain = 10, #TODO: get from metadata
            coverage_y_grain = 10,  #TODO: get from metadata
            coverage_t_grain = 1
        )


        ## TODO: account for different bands

        level = self.createNewLevel("initial_level")
        #
        # Note there are some special treatments since this is the loading procedure (hen-egg problem). Do not 
        # do this somewhere else. Basically it avoids that the atoms have to be created multiple times. The numpy array 
        # should not be further exposed to the user as the theory says it is actually a field of atoms.
        #
        array = ((area_of_interest.nir - area_of_interest.red) / (area_of_interest.nir + area_of_interest.red))

        atoms = self.extractAtoms(array, domain_shape)

        self.setCoverageView(level["depth"], array, atoms)
       
        self.createObjectViewFromCoverage(level["depth"], create_level=False, algorithm="pixelwise")

        #
        # Set default neighbourhood concept
        #
        Datacube.neighbourhood = Neighbourhood("4-connected")

    #
    # Calculate the Index position for the atom
    # TODO: Add time dimension
    #
    def getIndexPosition(self, coordinate):
        """
        Calculates the index position for the atom tuple.

        :param coordinate: tuple; coordinates
        :return: int
        """                
        x_coord = coordinate[0]
        y_coord = coordinate[1]
        t_coord = coordinate[2]

        if x_coord > self.i_xmax or x_coord < self.i_xmin:
            print ("invalid x coordinate: " + str(x_coord) + " because it is smaller than " + str(self.i_xmin) + " or larger than " + str(self.i_xmax))
            return None
        
        if y_coord > self.i_ymax or y_coord < self.i_ymin:
            print ("invalid y coordinate: " + str(y_coord))
            return None

        if t_coord > self.i_tmax or t_coord < self.i_tmin:
            print("invalid t coordinate")
            return None

        i_xgrain = (self.i_xmax - self.i_xmin) / self.i_number_of_cells
        i_ygrain = (self.i_ymax - self.i_ymin) / self.i_number_of_cells
        i_tgrain = (self.i_tmax - self.i_tmin) / self.i_number_of_cells

        x_hilbert = max(math.floor((x_coord - self.i_xmin) / i_xgrain), self.i_number_of_cells)
        y_hilbert = max(math.floor((y_coord - self.i_ymin) / i_ygrain), self.i_number_of_cells)
        t_hilbert = max(math.floor((t_coord - self.i_tmin) / i_tgrain), self.i_number_of_cells)

        return self.sfc.distance_from_coordinates([x_hilbert, y_hilbert, t_hilbert])


    def getIndex(self):
        return self.sfc_index

        
    def extractAtoms(self, target, domain_shape):

        atoms = []

        metadata = self.getDatasetMetadata()

        if type(target) == xarray.DataArray:

            if Datacube.config["show_progress"]:
                counter = 0
                size = metadata["coverage_x_size"] * metadata["coverage_y_size"]
                rhino_tools.progressBar(counter, size, prefix = 'Extract atoms:                     ', suffix = 'Complete', length = 50)
        

            for x_index in range(0, metadata["coverage_x_size"]):
                for y_index in range(0, metadata["coverage_y_size"]):

                    t_index = 0

                    coords_y = target.coords["x"].item(x_index)
                    coords_x = target.coords["y"].item(y_index)
                    coords_t = t_index

                    if Datacube.config["sfc_index"]:
                        h = self.getIndexPosition((coords_x, coords_y, coords_t))
                    else:
                        h = -1

                    if shapely_geometry.Point(coords_y, coords_x).within(domain_shape): #TODO: This will be repeated later, might be redundant, but I don't want to break the code

                        value = target.item((t_index, y_index, x_index))

                        #
                        # Create atom
                        #
                        a = Atom(
                                    {"lat": coords_x, "lon": coords_y},
                                    {"property": Datacube.atom_property, "value": value},
                                    {"x":x_index,"y":y_index},
                                    {"h":h}
                                )

                        #
                        # Append atom to list
                        #
                        atoms.append(a)

                        #
                        # Add to index
                        #
                        if Datacube.config["sfc_index"]:
                            self.sfc_index[h].append(a)

                    else:
                        target[t_index][y_index][x_index] = numpy.nan

                    if Datacube.config["show_progress"]: 
                        counter += 1                          
                        rhino_tools.progressBar(counter, size, prefix = 'Extract atoms:                     ', suffix = 'Complete', length = 50)

        else:
            pass
        
        return atoms


    def createNewLevel(self, name ="no name available", description = "no description available"):
        """
        Creates a new level, which is pre-condition for creating any object or coverage view.

        :param name: string; name of the level (optional)
        :param description: string; decription of the level (optional)
        """        
        depth = len(self.__levels)
        self.__levels.append({
            "depth": depth,
            "name": name,
            "description": description 
        })

        self.__view["coverage"].append(Coverage)
        self.__view["objects"].append([])

        return self.__levels[depth]


    def getLevel(self, depth):
        """
        Returns the level by the depth

        :param depth: integer
        :return: dict
        """
        for level in self.__levels:
            if level["depth"] == depth:
                return level
        return None


    def createObjectViewFromCoverage(self, from_level, create_level = True, algorithm = "pixelwise", parameters = None):

        """
        Creates the coverage view in the rhino datacube.

        :param from_level: int, level from which the coverage should be taken
        :param create_level: ...
        :param algorithm: ...
        :param parameters: ...
        """        
        if algorithm not in ["pixelwise", "aggregation"]:
            raise Exception("Algorithm " + str(algorithm) + " is unknown")
        
        coverage = self.getCoverageView(from_level)

        if algorithm == "pixelwise":

            atom_coverage = coverage.getCoverage()
            object_list = []

            #
            # Iterate through the atom coverage and create an object for each atom
            #
            for row in atom_coverage:
                for cell in row:

                    #
                    # If the cell is empty, continue and do not create an object
                    #
                    if cell is None:
                        continue

                    #
                    # Create an object and assign the atom in the cell to it
                    #
                    o = Object()
                    o.create([cell])
                    object_list.append(o)
            
            self.setObjectView(from_level, object_list)
        
        if algorithm == "aggregation":
            
            condition = Datacube.atom_property.getSimilarityFunction()

            #
            # STEP 1: Select all objects, which match the condition. If no objects
            # match the condition, return here.
            #
            # This function generates a new level and views
            #
            objects, coverage, level = self.selectByCondition(from_level, condition, create_level=create_level)

            if len(objects) == 0:
                print("No objects found") #TODO: Raise exception here
                return

            #
            # STEP 2: Aggregate and update views
            #

            #
            # Aggregate/cluster using the Link class, it updates the levels accordingly ##TODO: Not sure whether this is good to have it hidden implicitly there
            #
            Link("aggregation", self.__view, level["depth"])

            #
            # return information about the new level
            #
            return self.getObjectView(level["depth"]), self.getCoverageView(level["depth"]), level


    def setObjectView(self, level, objects):
        """
        Sets the object view of a specific level. Will override existing object view.

        :param level: int, target level
        :param objects: list of objects
        """                
        self.__view["objects"][level] = objects


    def getObjectView(self, from_level):
        """
        Returns the object view of a specific level.

        :param from_level: int, target level
        :return: dict
        """                
        return self.__view["objects"][from_level]


    def createCoverageViewFromObjectView(self, from_level, create_level=False, algorithm = None, parameters = None):
        """
        Creates the coverage view in the rhino datacube.
        This can also be seen as masking, e.g., after selection of objects.

        :param from_level: level from which the objects should be taken
        :param to_level: target level
        :param algorithm: None (currently unused)
        :param parameters: None (currently unused)
        """
        metadata = self.getDatasetMetadata()
        old_stupid_array = self.__view["coverage"][0].getStupidArray()
        x_size = metadata["coverage_x_size"]
        y_size = metadata["coverage_y_size"]

        new_stupid_array = old_stupid_array.copy(deep=True)
        new_stupid_array = new_stupid_array * 0 #TODO initialise with nan
        
        atoms = []
        for o in self.__view["objects"][from_level]:
            for a in o.getAtoms():
                x,y = a.getIndex()
                value = a.getObservationValue()
                new_stupid_array[0][y][x] = value
                atoms.append(a)

        coverage = Coverage(x_size, y_size)
        coverage.create(atoms=atoms, array=new_stupid_array)

        self.__view["coverage"][from_level] = coverage


    def setCoverageView(self, level, array, atoms):
        """
        Sets the coverage of a specific level. Will override existing coverage.

        :param level: int, target level
        :param array: array (i.e. stupid array)
        :param atoms: atoms
        :return: dict
        """        

        metadata = self.getDatasetMetadata()
        coverage = Coverage(metadata["coverage_x_size"], metadata["coverage_y_size"])
        coverage.create(array=array, atoms=atoms)
        self.__view["coverage"][level] = coverage


    def getCoverageView(self, from_level):
        """
        Returns the coverage view of a specific level.

        :param from_level: int, target level
        :return: dict
        """        
        return self.__view["coverage"][from_level]
    

    def getDatasetMetadata(self):
        """
        Returns the metadata of the selected dataset.

        :return: dict
        """
        return Datacube.dataset_metadata


    def setDatasetMetadata(self, coverage_x_size, coverage_y_size, coverage_x_grain, coverage_y_grain):
        
        Datacube.dataset_metadata = {
            "coverage_x_size": coverage_x_size,
            "coverage_y_size": coverage_y_size,
            "coverage_x_grain": coverage_x_grain,
            "coverage_y_grain": coverage_y_grain
        }

    


    ## TODO: This is currently not possible as switching to another neighbourhood
    #        would invalidate existing object links.
    #def setNeighbourhoodConcept(self, concept):
    #    neighbourhood = Neighbourhood(concept)
    

    def selectByCondition(self, from_level, condition, create_level = False):
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

        if Datacube.config["show_progress"]:
            counter = 0
            size = len(candidates)
            rhino_tools.progressBar(counter, size, prefix = 'Selecting objects with condition:  ', suffix = 'Complete', length = 50)

        for obj in candidates:

            if locals()[condition["operator"]](obj.getAttributeValue(condition["property"]), condition["value"]):
                objects.append(obj)

            if Datacube.config["show_progress"]:
                counter += 1
                rhino_tools.progressBar(counter, size, prefix = 'Selecting objects with condition:  ', suffix = 'Complete', length = 50)
        
        if create_level == True:

            level = self.createNewLevel()
            self.__view["objects"][level["depth"]] = objects
            self.createCoverageViewFromObjectView(level["depth"])
            coverage = self.__view["coverage"][level["depth"]]

        else:

            level = None
            coverage = None
    
        return objects, coverage, level


    def execute(self, procedure, parameters=None, level = 0):
        """
        This is a wrapper to execute a procedure either on the object-based or on the coverage-based view.
        It creates a new level by default

        :param procedure: procedure to be executed
        :param parameter: optional parameter for the procedure
        """
        existing_procedures = ["aggregation", "selection"]

        if procedure not in existing_procedures:
            raise Exception("Choose one of " + ",".join(existing_procedures))
        
        if procedure == "aggregation":
            return self.createObjectViewFromCoverage(level, algorithm = "aggregation", parameters = None, create_level = True)
        
        if procedure == "selection":
            return self.selectByCondition(level, parameters, create_level = True)

    def executeQuery(self, command):
        """
        Executes a SQL-like query command against the object view of the rhino data cube.

        :param command: string
        :return: ?
        """
        pass ##TODO: Implement


class Atom:

    def __init__(self, coordinates, observation, index, h):

        ## TODO: Check value integrity and trow exceptions
   #     rhino_tools.checkCoordinates(coordinates["lat"], coordinates["lon"])

        self.__tuple = {
            "id": id(self),
            "lat": coordinates["lat"],
            "lon": coordinates["lon"],
            "geometry": shapely_geometry.Point(coordinates["lon"],coordinates["lat"]),
            "property": observation["property"],
            "value": observation["value"],
            "x_index": index["x"],
            "y_index": index["y"],
            "h": h
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

    def _print(self):
        print(self.__tuple)

class Coverage:

    def __init__(self, x_size, y_size):
        self.id = id(self)
        self.coverage = []
        self.__array = None

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
            self.__array = array


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
        return self.__array


    def getSize(self):
        """
        Returns the size in x and y direction of the coverage.

        :return: tuple (x_size, y_size)
        """                 
        return (self.x_size, self.y_size)


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
        self.__calculateGeoAttributes()

    def grow(self, atoms):
        """
        Grows the object with the given atomes

        :param atoms: list of atoms
        """
        self.atoms.extend(atoms)
        self.__calculateGeoAttributes()
        self.__attributes["derived"]["geometry"] = None

    def shrink(self, atoms):
        """
        Shrinks the object by the given atomes

        :param atoms: list of atoms
        """
        self.__attributes["derived"]["geometry"] = None

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
            # TODO: Search for better algorithm to delineate border
            #
            buf = Datacube.neighbourhood.getMaxDistance()
            geometry = shapely_operations.cascaded_union(multipoints.buffer(buf))
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

            if attribute == "area":
                geom = self.getGeometry()
                return geom.area

        else:
                
            matches = re.match( r'^((avg)|(min)|(max)|(std))\(([a-z]+)\)$', attribute, re.I)
            aggregate_function = matches.group(1)
            property_type = matches.group(6)

            if aggregate_function is not None and property_type in self.__attributes["modelled"]:

                atomProperties = numpy.array([])
                aggregate_value = None

                for a in self.atoms:
                    atomProperties = numpy.append(atomProperties, a.getObservationValue()) ## TODO: return based on the property arg
                
                if len(atomProperties) == 1:
                    return atomProperties[0]

                if aggregate_function == "avg":
                    aggregate_value = numpy.average(atomProperties)
                if aggregate_function == "min":
                    aggregate_value = numpy.min(atomProperties)
                if aggregate_function == "max":
                    aggregate_value = numpy.max(atomProperties)
                if aggregate_function == "std":
                    aggregate_value = numpy.std(atomProperties)

                return aggregate_value

            else:
                pass

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


    def _print(self):
        """
        Prints the object.
        """                
        print("{ id: " + str(self.id) + ", atoms: ") 
        for a in self.atoms:
            a._print()
        print(",links: ")
        for l in self.links:
            print(l)
        print("}")

class Link:

    def __init__(self, link_type, target, level, mutual = True):

        self.id = id(self)
        self.attributes = copy.deepcopy(Datacube.link_model)

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

        array = coverage.getStupidArray()
        atom_coverage = coverage.getCoverage(None)

        x_size, y_size = coverage.getSize()

        #
        # Get regions using the stupid array
        #
        from skimage.measure import label as region_group

        grouped = region_group((array/array).fillna(0), background=0)
        #
        # Create a nested list where the atoms will be stored
        #
        atom_lists = []
        for i in range(grouped.max()):
            atom_lists.append([])

        #
        # Aggregate the atoms to the atom_list
        #
        if Datacube.config["show_progress"]:
            counter = 0
            rhino_tools.progressBar(counter, x_size*y_size, prefix = 'Object linking (aggregate):        ', suffix = 'Complete', length = 50)

        for x in range(x_size):
             for y in range(y_size):

                #
                # Get the object number of this location
                #
                object_number = grouped[0][y][x]

                #
                # We treat objects with number 0 as background
                #
                if object_number != 0:

                    #
                    # Get the atom of this location
                    #
                    atom = atom_coverage[x][y]
                    if atom is None:
                        print("Fail: No atom found, but there should be one") #TODO: Make exception
                    else:
 #                       print(atom)
 #                       #print(atom._print())
#
                        atom_lists[object_number-1].append(atom)

                if Datacube.config["show_progress"]:
                    counter += 1
                    rhino_tools.progressBar(counter, x_size*y_size, prefix = 'Object linking (aggregate):        ', suffix = 'Complete', length = 50)

        #
        # Iterate through the atom list and create objects
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


class Property:

    def __init__(self, property_type):
        self.property_type = property_type
        self.aggregation_function = "any"


    def setSimilarityFunction(self, similarity):
        self.similarity = similarity
    

    def getSimilarityFunction(self):
        return self.similarity


    def setAggregationFunction(self, function):
        if function not in ["any","all"]:
            raise Exception("Aggregation function must be one of any, all")
        self.aggregation_function = function


    def getPropertyType(self):
        return self.property_type


class Neighbourhood:

    def __init__(self, concept = "4-connected"):
        
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