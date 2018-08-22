import rhino
import datetime
import png

f = open("demo-run.txt","w")
f.write("Demo-data execution start: ")
f.write(str(datetime.datetime.now()) + "\n")

f.write("1: Load demo dataset ...")
dc = rhino.Datacube(show_progress = True)
dc.load('demodata/demolayer_lowres.tif', None)
f.write("done\n")


f.write("2: Aggregation based on colour information ...")
objects, coverage, level = dc.execute("aggregation", parameters = {"key":"class","value":1, "operator": "eq"}, level = 0)
f.write("done\n")
png.from_array(coverage.getStupidArray().tolist(), 'L').save("result-level-"+str(level["depth"])+".png")

f.write("3: Selection based on shape information ....")
objects, coverage, level = dc.execute("selection", parameters = {"key": "compactness", "value": 0.5, "operator": "gt"}, level= level["depth"])
f.write("done\n")
png.from_array(coverage.getStupidArray().tolist(), 'L').save("result-level-"+str(level["depth"])+".png")
f.write("Demo-data execution finished.")

f.close()




