from rhino import Datacube
from rhino import Atom

import rhino_helper

from shapely import geometry as shapely_geometry
import unittest

class TestRhinoHelper(unittest.TestCase):

    def test_checkcoordinates(self):

        lat = 0.0
        lon = 0.0
        self.assertTrue(rhino_helper.checkCoordinates(0,0))

        lat = 100.0
        lon = 0.0
        with self.assertRaises(rhino_helper.InvalidCoordinatesError) as context:
            rhino_helper.checkCoordinates(lat, lon)
        self.assertEqual("latitude " + str(lat) + " not allowed", str(context.exception))

        lat = -100.0
        lon = 0.0
        with self.assertRaises(rhino_helper.InvalidCoordinatesError) as context:
            rhino_helper.checkCoordinates(lat, lon)
        self.assertTrue("latitude " + str(lat) + " not allowed" == str(context.exception))

        lat = 0.0
        lon = 200.0
        with self.assertRaises(rhino_helper.InvalidCoordinatesError) as context:
            rhino_helper.checkCoordinates(lat, lon)
        self.assertTrue("longitude " + str(lon) + " not allowed" == str(context.exception))

        lat = 0.0
        lon = -200.0
        with self.assertRaises(rhino_helper.InvalidCoordinatesError) as context:
            rhino_helper.checkCoordinates(lat, lon)
        self.assertTrue("longitude " + str(lon) + " not allowed" == str(context.exception))

        lat = 100.0
        lon = 200.0
        with self.assertRaises(rhino_helper.InvalidCoordinatesError) as context:
            rhino_helper.checkCoordinates(lat, lon)
        self.assertTrue("latitude " + str(lat) + " not allowed" == str(context.exception))


# class TestDatacubes(unittest.TestCase):

#     def setUp(self):
#         self.dc = Datacube()
#         self.dc.load('demodata/demolayer_lowres.tif', None)

#     def test_load(self):
#         assert self.dc.dataset is not None
    
#     def test_level(self):
#         # initial level
#         level = {
#             "level": 0,
#             "name": "initial_level",
#             "description": "no description available"
#             }

#         self.assertDictEqual(self.dc.getLevel(0), level)

#         new_level = {
#             "level": 1,
#             "name": "New test level",
#             "description": "This is the description for the new test level"
#         }
#         self.dc.createNewLevel(new_level["name"],new_level["description"])

#         self.assertDictEqual(self.dc.getLevel(1), new_level)

#     def test_aggregateandselectobjects(self):
#         level = self.dc.aggregate(0, {"key":"class","value":1, "operator": "eq"})
#         objects, coverage, level = self.dc.selectObjectsByCondition(level["level"],{"key": "compactness", "value": 0.5, "operator": "gt"}, create_level=True)
#         self.assertEqual(len(objects),1)
#         self.assertEqual(level,2)

#     def test_getdatasetmetadata(self):
#         metadata = self.dc.getDatasetMetadata()
#         self.assertDictEqual(metadata, {'coverage_x_size': 100, 'coverage_y_grain': -0.020677661659999985, 'coverage_y_size': 100, 'coverage_x_grain': 0.02314663617999997})

class TestAtoms(unittest.TestCase):
    
    def setUp(self):
        self.latitude = 0
        self.longitude = 0
        self.coordinates = {"lat": self.latitude, "lon": self.longitude}
        self.observation = {"property": "class", "value":0}
        self.index = {"x":0,"y":0}

        self.a = Atom(self.coordinates, self.observation, self.index)

    def test_getid(self):
        self.assertGreaterEqual(self.a.getID(),0)

    def test_getcoordinates(self):
        _tested_coordinates = self.a.getCoordinates()

        self.assertEqual(_tested_coordinates[0],self.latitude)
        self.assertEqual(_tested_coordinates[1],self.longitude)

    def test_getlatititude(self):
        self.assertEqual(self.a.getLatitude(),self.latitude)

    def test_getlongitude(self):
        self.assertEqual(self.a.getLongitude(),self.longitude)

    def test_getObservationValue(self):
        _tested_observation = self.a.getObservationValue()
        
        self.assertEqual(_tested_observation, self.observation["value"])

    def test_getObservationProperty(self):
        _tested_observation = self.a.getObservationProperty()
        
        self.assertEqual(_tested_observation, self.observation["property"])

    def test_getObservation(self):
        _tested_observation = self.a.getObservation()
        self.assertDictEqual(self.observation, _tested_observation)

    def test_getgeometry(self):
        _tested_geometry = self.a.getGeometry()
        self.assertEqual(shapely_geometry.Point(self.latitude, self.longitude),_tested_geometry)

    def test_getindex(self):
        _tested_index = self.a.getIndex()

        self.assertEqual(self.index["x"],_tested_index[0])
        self.assertEqual(self.index["y"],_tested_index[1])

if __name__ == '__main__':
    unittest.main()