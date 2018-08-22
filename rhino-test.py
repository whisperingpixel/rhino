from rhino import Datacube
from rhino import Atom
from rhino import Object
from rhino import Coverage

import rhino_helper

from shapely import geometry as shapely_geometry
import gdal
import numpy

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


class TestDatacubes(unittest.TestCase):

    def setUp(self):
        self.dc = Datacube()
        self.dc.load('demodata/demolayer_lowres.tif', None)

    def test_load(self):
        assert self.dc.dataset is not None

    def test_level(self):
        #initial level
        level = {
            "depth": 0,
            "name": "initial_level",
            "description": "no description available"
            }

        self.assertDictEqual(self.dc.getLevel(0), level)

        new_level = {
            "depth": 1,
            "name": "New test level",
            "description": "This is the description for the new test level"
        }
        self.dc.createNewLevel(new_level["name"],new_level["description"])

        self.assertDictEqual(self.dc.getLevel(1), new_level)

    def test_aggregateandselectobjects(self):
        objects, coverage, level = self.dc.createObjectViewFromCoverage(0, algorithm = "aggregation", parameters = {"key":"class","value":1, "operator": "eq"}, create_level=True, )
        objects, coverage, level = self.dc.selectObjectsByCondition(level["depth"],{"key": "compactness", "value": 0.5, "operator": "gt"}, create_level=True)
        self.assertEqual(len(objects),1)
        self.assertEqual(level["depth"],2)

    def test_getdatasetmetadata(self):
        metadata = self.dc.getDatasetMetadata()
        self.assertDictEqual(metadata, {'coverage_x_size': 100, 'coverage_y_grain': -0.020677661659999985, 'coverage_y_size': 100, 'coverage_x_grain': 0.02314663617999997})

    def test_getcoverageview(self):

        testdataset = gdal.Open("demodata/demolayer_lowres.tif")
        band = testdataset.GetRasterBand(1)
        numpy_array = numpy.array(band.ReadAsArray())

        coverage = self.dc.getCoverageView(0)
        coverage_array = coverage.getStupidArray()

        for (x_index, y_index), value in numpy.ndenumerate(numpy_array):
            self.assertEqual(coverage_array[x_index][y_index], value)

    def test_demodata_workflow_coverage(self):
        #
        # STEP 1: Import data. The coverage should be exactly the same as the input.
        #
        testdataset = gdal.Open("demodata/demolayer_lowres_result_level_0.tif")
        band = testdataset.GetRasterBand(1)
        numpy_array = numpy.array(band.ReadAsArray())

        coverage = self.dc.getCoverageView(0)
        coverage_array = coverage.getStupidArray()

        for (x_index, y_index), value in numpy.ndenumerate(numpy_array):
            self.assertEqual(coverage_array[x_index][y_index], value)

        #
        # STEP 2: Aggregate. The coverage should cover only pixels with class value==1.
        #
        objects, coverage, level = self.dc.createObjectViewFromCoverage(0, algorithm = "aggregation", parameters = {"key":"class","value":1, "operator": "eq"}, create_level=True)

        testdataset = gdal.Open("demodata/demolayer_lowres_result_level_1.tif")
        band = testdataset.GetRasterBand(1)
        numpy_array = numpy.array(band.ReadAsArray())

        coverage = self.dc.getCoverageView(level["depth"])
        coverage_array = coverage.getStupidArray()

        for (x_index, y_index), value in numpy.ndenumerate(numpy_array):
            self.assertEqual(coverage_array[x_index][y_index], value)

        #
        # STEP 3: Select: The coverage should contain only pixels with class value == 1 and the object's compactness value >= 0.5.
        #
        objects, coverage, level = self.dc.selectObjectsByCondition(level["depth"],{"key": "compactness", "value": 0.5, "operator": "gt"}, create_level=True)

        testdataset = gdal.Open("demodata/demolayer_lowres_result_level_2.tif")
        band = testdataset.GetRasterBand(1)
        numpy_array = numpy.array(band.ReadAsArray())

        coverage = self.dc.getCoverageView(level["depth"])
        coverage_array = coverage.getStupidArray()

        for (x_index, y_index), value in numpy.ndenumerate(numpy_array):
            self.assertEqual(coverage_array[x_index][y_index], value)


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


class TestObjects(unittest.TestCase):

    def setUp(self):

        self.dc = Datacube()
        self.dc.load('demodata/demolayer_lowres.tif', None)
        objects, coverage, level = self.dc.createObjectViewFromCoverage(0,  algorithm = "aggregation", parameters = {"key":"class","value":1, "operator": "eq"},create_level=True)
        objects, coverage, level = self.dc.selectObjectsByCondition(level["depth"],{"key": "compactness", "value": 0.5, "operator": "gt"}, create_level=True)
        self.object2 = objects[0]

        self.atoms = []
        self.atoms.append(
            Atom({"lat": 0, "lon": 0},{"property": "class", "value":0},{"x":0,"y":0}))
        self.atoms.append(
            Atom({"lat": 0, "lon": 10},{"property": "class", "value":0},{"x":0,"y":10}))
        self.atoms.append(
            Atom({"lat": 10, "lon": 10},{"property": "class", "value":0},{"x":10,"y":10}))
        self.atoms.append(
            Atom({"lat": 10, "lon": 0},{"property": "class", "value":0},{"x":10,"y":0}))

        self.object1 = Object()
        self.object1.create(self.atoms)

    def test_create(self):

        _tested_atoms = self.object1.atoms

        self.assertEqual(len(self.atoms), len(_tested_atoms))

        for i,a in enumerate(_tested_atoms):
            coords1 = a.getCoordinates()
            coords2 = self.atoms[i].getCoordinates()
            self.assertEqual(coords1, coords2)


    def test_getid(self):
        self.assertGreaterEqual(self.object1.getID(),0)


    def test_getatoms(self):

        _tested_atoms = self.object1.getAtoms()

        self.assertEqual(len(self.atoms), len(_tested_atoms))

        for i,a in enumerate(_tested_atoms):
            coords1 = a.getCoordinates()
            coords2 = self.atoms[i].getCoordinates()
            self.assertEqual(coords1, coords2)


    def test_getnumberofatoms(self):

        _tested_atoms_length = self.object1.getNumberOfAtoms()

        self.assertEqual(len(self.atoms), _tested_atoms_length)


    def test_getatomcoordinates(self):
        self.assertEqual(
            self.object1.getAtomCoordinates(),
            [(0, 0), (0, 10), (10, 10), (10, 0)]
        )


    def test_getboundingbox(self):
        self.assertEqual(
            self.object1.getBoundingBox(),
            shapely_geometry.Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)])
        )


    def test_getgeometry(self):
        pass


    def test_getattributevalue(self):
        self.assertAlmostEqual(self.object2.getAttributeValue("compactness"),0.7219841742616753)

if __name__ == '__main__':
    unittest.main()