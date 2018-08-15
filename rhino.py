import gdal
import numpy
import shapely
import math

class Datacube:

    dataset = None

    object_model = {
        "derived": ["number_of_atoms"],
        "modelled": ["class"]
    }

    link_model = {
        "derived": ["link_type"],
        "modelled": []
    }

    neighbourhood = None

    def __init__(self):
        self.__datacube_coverage = []
        self.__datacube_objects = []


    def load(self, boundingbox):
        """
        Loads the rhino datacube.

        :param boundingbox: [x_min, x_max, y_min, y_max]
        """
        ## TODO: Load datacube from here. Be advised that this is a workaraound without a datacube

        Datacube.dataset = gdal.Open('demodata/demolayer.tif')
        
        ## TODO: account for different bands
        band = Datacube.dataset.GetRasterBand(1)
        self.__datacube_coverage = numpy.array(band.ReadAsArray())

        #
        # Set default neighbourhood concept
        #
        Datacube.neighbourhood = Neighbourhood("4-connected")


    def createObjectView(self, layer):
        """
        Creates the object view in the rhino datacube.

        :param layer: int
        """
        (upper_left_x, x_size, x_rotation, upper_left_y, y_rotation, y_size) = Datacube.dataset.GetGeoTransform()

        #
        # Iterate through all pixels and create an atom from it. Then, create an object for each atom.
        # Initial objects only consist of one single atom
        #
        # Stolen from here: https://stackoverflow.com/questions/6967463/iterating-over-a-numpy-array/6967491#6967491
        for (x_index, y_index), value in numpy.ndenumerate(self.__datacube_coverage):
            
            #
            # Extract coordinates, properties and values from atoms
            # Stolen from here: https://gis.stackexchange.com/questions/42790/gdal-and-python-how-to-get-coordinates-for-all-cells-having-a-specific-value/42846#42846
            x_coord = x_index * x_size + upper_left_x + (x_size / 2) #add half the cell size
            y_coord = y_index * y_size + upper_left_y + (y_size / 2) #to centre the point
            # TODO: add time
            property = "class"
            value = int(self.__datacube_coverage[x_index, y_index]) ## TODO: remove int

            #
            # Create an object and fill it with the atom
            #
            o = Object()
            o.create([Atom({"lat": x_coord, "lon": y_coord},{"property": property, "value": value})])
            self.__datacube_objects.append(o)


    def createCoverageView(self, objects):
        """
        Creates the coverage view in the rhino datacube.

        :param objects: list of objects
        """        
        pass

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
    def selectObjectsByCondition(self, condition):
        """
        Select the objects in the object view of the rhino datacube, which match
        a certain condition.

        :param condition: ?
        :return: list of objects
        """        
        objects = []
        candidates = self.__datacube_objects
        for obj in candidates:
            if obj.getAttributeValue("class") == 1:
                objects.append(obj)
        return objects


    def aggregateObjects(self, condition):
        """
        Aggregates objects based on a certain condition.

        :param condition: ?
        :return: ?
        """      
        #
        # Get all objects, which match the condition. If no objects
        # match the condition, return here.
        #
        objects = self.selectObjectsByCondition(condition)
        if len(objects) == 0:
            return

        #
        # Aggregate/cluster the objects using the object link
        # class.
        #
        l = Link("aggregate", objects, True)
        for o in l.getObjects():
            o.print()


    def executeQuery(self, command):
        """
        Executes a SQL-like query command against the object view of the rhino data cube.

        :param command: string
        :return: ?
        """
        pass ##TODO: Implement


class Atom:

    def __init__(self, coordinates, observation):

        ## TODO: Check value integrity
        self.__tuple = {
            "id": id(self),
            "lat": coordinates["lat"],
            "lon": coordinates["lon"],
            "property": observation["property"],
            "value": observation["value"]
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
    
    def print(self):
        print(self.__tuple)


class Object:

    def __init__(self):
        self.id = id(self)
        self.atoms = []
        self.links = []
        self.attributes = {"derived":{}, "modelled":["class"]} ## TODO: This should be derived from the object model

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
        self.atoms.extend(atoms)
#        self.__calculateAlphaShape(3)
        self.__calculateGeoAttributes()

    def shrink(self, atoms):
        pass

    def __calculateAlphaShape(self, neighbors):
        
        from rhino_geo import calcualateAlphaShape
        polygon = calcualateAlphaShape(self.getAtomCoordinates(), neighbors)
        if polygon:
            self.attributes["derived"]["outline"] = polygon
        ## TODO: Handle case where not enough atoms are available for that polygon

    def __calculateGeoAttributes(self):
        self.attributes["derived"]["number_of_atoms"] = len(self.atoms)

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
        if attribute in self.attributes["derived"]:
            return self.attributes["derived"][attribute]
        else:
            if attribute in self.attributes["modelled"]:
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
        if attribute in self.attributes["modelled"]:
            self.attributes["modelled"][attribute] = value
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

    def __init__(self, link_type, objects,  mutual):
        self.id = id(self)
        self.attributes = {"derived":{}}
        self.attributes["derived"]["link_type"] = link_type

        self.__new_objects = []

        if link_type == "aggregate":
            
            self.__new_objects = self.__aggregate(objects, Datacube.neighbourhood)

            # demo_array = self.__demo_aggregate(objects)
            # print("end clustering")
            # for atoms in demo_array:
            #     #atoms = []
            #     #for o in objects:
            #     #    atoms.extend(o.getAtoms())

            #     new_o = Object()
            #     new_o.create(atoms)

            #     self.__new_objects.append(new_o)

        else:
            if len(object != 2):
                raise Exception ("Exactly 2 objects are expected to be part of the link")

            target = objects[0]
            candidate = objects[1]

            target.linksWithObjects(link_type, candidate)
            if mutual == True:
                candidate.linksWithObjects(link_type, target)

    def __aggregate(self, objects, neighbourhood):
        """
        Aggregates the selected objects using the given concept of neighbourhood

        :param objects: List of objects to be aggregated
        :param neighbourhood: Neighbourhood object (i.e. inst. class)
        :return: list of (aggregated) objects
        """

        #
        # Position of the list, initial starting point will be 0
        #
        list_pos = 0

        #
        # As we are juggling with the object list, make a copy of it
        #
        aggregates = objects.copy()
        final_aggregates = []

        #
        # Get the first object
        #
        temp_obj = aggregates[0]
        #
        # Remove this object from the list, so that it will not be used for neighborhood estimation
        #
        del aggregates[0]

        #
        # Iterate through the list of aggregates
        #
        while len(aggregates) > 0:

            #
            # Set a flag whether we found something or not
            #
            aggregated = False

            #
            # Iterate through the other objects
            #

            for i,cand in enumerate(aggregates):
                #
                # If the object is neighbour
                #
                if neighbourhood.isNeighbour(temp_obj, cand):

                    #
                    # Grow the current object
                    #
                    temp_obj.grow(cand.getAtoms())

                    #
                    # delete the candidate from the list
                    #
                    del aggregates[i]

                    #
                    # Set the flag that a neighbour has been found
                    #
                    aggregated = True
            
            if aggregated == False:

                #
                # put the object back in the list if we are finished
                #
                final_aggregates.append(temp_obj)

                #
                # increase list position if we did not find neighbours
                #                
                list_pos += 1

                #
                # Get the next object
                #
                temp_obj = aggregates[0]                

        return final_aggregates

    def __demo_aggregate(self, objects):

        from sklearn.cluster import DBSCAN
        np_arr = []
        at_arr = []
        for o in objects:
            atoms = o.getAtoms()
            for a in atoms:
                at_arr.append(a)
                np_arr.append([a.getLatitude(), a.getLongitude()])

        X = numpy.array(np_arr)
        db = DBSCAN(eps=0.005).fit(X)

        atom_list = []
        for i,label in enumerate(db.labels_):
            try:
                atom_list[label].append(at_arr[i])
            except Exception:
                atom_list.append([at_arr[i]])
        
        print (len(atom_list))
        return atom_list

    def getObjects(self):
        """
        Returns the objects involved in the link (relationship/aggregation)

        :return: list of objects
        """                
        return self.__new_objects

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

    def isNeighbour(self, object, candidate):
        """
        Compares two objects and evaluates whether they are neighbours or not, based on the
        selected concept of neighbourhood.

        :param object: Target object
        :param candidate: 
        :return: boolean
        """                
        min_dist = float("inf")

        obj_atoms = object.getAtoms()
        cand_atoms = candidate.getAtoms()

        for a1 in obj_atoms:
            for a2 in cand_atoms:
                dist = math.hypot(a1.getLatitude() - a2.getLatitude(), a1.getLongitude() - a2.getLongitude())
                if dist < min_dist:
                    min_dist = dist

        ## TODO: Account for different grains in X and Y direction

        #
        # Get allowed distance
        #
        allowed_dist = max(Datacube.getDatasetMetadata(Datacube)["coverage_x_grain"], Datacube.getDatasetMetadata(Datacube)["coverage_y_grain"])

        #
        # add tolerance
        #
        allowed_dist += self.tolerance

        #
        # Check whether the detected distance is smaller than the allowed distance
        #
        if min_dist < allowed_dist: ##TODO: this is for testing
            return True
        else:
            return False
