from rhino import rhino

atom_property = rhino.Property("ndvi")
atom_property.setSimilarityFunction({"property":"min(ndvi)","value":0.2, "operator": "gt"})
dc = rhino.Datacube(show_progress = True, atom_property=atom_property)

dc.load(product = "s2a_sen2cor_granule_ingested_10", domain = {"x":(676388, 680313), "y":(2453906,2459794), "crs":'EPSG:32635'})


#dc.load(product = "s2a_sen2cor_granule_ingested_10", boundingbox = {"x":(676388, 677000), "y":(2453906,2454000), "crs":'EPSG:32635'})


objects, coverage, level = dc.execute("aggregation", parameters = None, level = 0)

objects[0].getAttributeValue("avg(ndvi)")