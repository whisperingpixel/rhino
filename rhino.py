import math

class Datacube:

    def __init__(self):
        self.datacube_coverage = []
        self.datacube_object = []

        self.createCoverageView()
        self.createObjectView()
    
    def createCoverageView(self):
        ## TODO: Load datacube from here
        self.datacube_coverage = None

    def createObjectView(self):
        
        for cell in self.datacube_coverage:
            o = Object()
            ##TODO: Extract values from cell

            lat = None
            lon = None
            property = None
            value = None

            
            o.create([Atom({"lat": lat, "lon": lon},{"property": property, "value": value})])
            self.datacube_object.append(o)

    def executeQuery(self, command):
        pass ##TODO: Implement

    #
    # These are user-friendly wrapper-functions
    #
    def aggregateObjects(self, condition):
        pass
    
    def createContract(self, condition):
        pass
        
class Atom:

    def __init__(self, coordinates, observation):

        ## TODO: Check value integrity
        self.tuple = {
            "id": id(self),
            "lat": coordinates["lat"],
            "lon": coordinates["lon"],
            "property": observation["property"],
            "value": observation["value"]
        }

    def euclideanDistanceTo(self,candidate):
        return math.hypot(self.tuple["lat"] - candidate.tuple["lat"], self.tuple["lon"] - candidate.tuple["lon"])

    def getObservationProperty(self):
        return self.tuple["property"]

    def getObservationValue(self):
        return self.tuple["value"]

    def getObservation(self):
        return {"property": self.tuple["property"], "value": self.tuple["value"]}

    def print(self):
        print(self.tuple)


class Object:

    def __init__(self):
        self.id = id(self)
        self.atoms = []
        self.links = []

    def create(self, atoms):
        self.atoms = atoms
    
    def getAtoms(self):
        return self.atoms

    def linksWithObjects(self, link_type, candidate):

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

    def print(self):
        print("{ id: " + str(self.id) + ", atoms: ") 
        for a in self.atoms:
            a.print()
        print(",links: ")
        for l in self.links:
            print(l)
        print("}")


class Link:

    def __init__(self, link_type, objects, condition, mutual):
        self.id = id(self)
        self.type = link_type
        
        if link_type == "aggregate":
            self.new_objects = []

            atoms = []
            for o in objects:
                atoms.extend(o.getAtoms())

                if condition:
                    pass

            new_o = Object()
            new_o.create(atoms)
            self.new_objects.append(new_o)

        else:
            if len(object != 2):
                raise Exception ("Exactly 2 objects are expected to be part of the link")

            target = objects[0]
            candidate = objects[1]

            target.linksWithObjects(link_type, candidate)
            if mutual == True:
                candidate.linksWithObjects(link_type, target)

    def getObjects(self):
        return self.new_objects