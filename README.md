![](https://travis-ci.org/whisperingpixel/rhino.svg?branch=master) [![codecov](https://codecov.io/gh/whisperingpixel/rhino/branch/master/graph/badge.svg)](https://codecov.io/gh/whisperingpixel/rhino) [![Documentation Status](https://readthedocs.org/projects/rhino-doc/badge/?version=latest)](https://rhino-doc.readthedocs.io/en/latest/?badge=latest)

# RHINO

## What is RHINO?

RHINO is a translation layer, which allows switching between a coverage view and an object view within a data cube.

## What does RHINO mean?

The name RHINO does not mean anything (or everything, depending on your attitude and whether you are more optimistic or pessimistic). To be honest, I chose it because it is short and easy to remember and you can make a nice logo from it:

![rhino_logo](https://raw.githubusercontent.com/whisperingpixel/rhino/master/res/logo.png)



## Why do I want to have this?

Geographers (and people from other domains, too) traditionally like to separate their view of the world either into a *coverage*-based or *object*-based representation (sometimes this is referred to as raster and vector, although it is not quite the same). This is documented by a number of software and approaches, which support usually only one of them. However, there are several application areas, where an analyst has to switch between both worlds seamlessly. RHINO solves this gap by building a bridge between the two worlds. It does not matter whether you want to perform object-based or coverage-based analysis, in RHINO you can conduct both analyses and can even combine them into an integrated workflow. 

Now, for the sake of completeness, it is true that modern GIS or remote-sensing software (ArcGIS, eCognition) allows switching between the two worlds to some degree. However, RHINO is especially designed to be built on top of data cubes and implements some of the more recent considerations about the nature of geospatial information coming from GIScience theory.

## How does it work?

When you access a data cube such as the Open Data Cube, RHINO internally creates an object-based view and a coverage-based view on the first level and keeps them synchronised. Each step of an analysis might go to either of the views, depending on the intention of the analyst, what is more efficient, or what is possible at all in a specific representation. During the analysis it is not necessary to stick to one representation, e.g, the first analysis step might be in the coverage-based view (such as: "aggregate all water cells"), the second step might proceed in the object-based view (such as: "select water objects having a certain size and shape"), and the third step might be conducted in the coverage-based view again (such as: "calculate the average temperature using the thermal band"). For each step, a new level will be generated, which allows you to go back and forth any time. The beauty is that you don't have to care about creating and managing the views or keeping them synchronised. This is all done by RHINO, you just need to conduct the analysis. 

The most basic element in RHINO is the **atom** (Goodchild, Yuan & Cova, 2017). An atom relates a location on the Earth to a property and its value. For example at location [lat] / [lon] the property *land cover* might be *water*. This is what an atom can tell and it is the only thing it can do content-wise (it may have additional qualities such as grain, accuracy etc). It is mainly driven by an observation point-of-view, where the observations are usually point-based (e.g. with a remote satellite sensor or a in-situ ocean buoy). Both of the two representations of the world make use of the same atoms, but it depends on your preferred representation of the world what you "do" with them.

You can put all of the the atoms into a regular or irregular grid.  You chose to view the atoms in a **coverage**-based representation. (We treat the grid as a synonym of a coverage for now, although strictly it is not correct).

On the other hand, you can also put all of the surrounding atoms of the aforementioned water atom, which have the same property value, into a common "container" and draw a border around them. Now you chose to have an **object**-based representation. An object has properties such as size, shape, etc and fulfils certain criteria regarding geometry and homogeneity (Blaschke, 2010, Lang et al. 2014).

In the physical world, objects are connected to each other in a explicit or implicit network, they maintain a relationship (association) between them (Kuhn, 2012). In RHINO the relationships or association between objects are called **link**. The reason is that linking objects might generate a new object (e.g. by generalisation of surface water bodies) as well as maintain a relationship between them (e.g. a contract between a cloud and cloud shadow).

Determining which objects are linked (e.g. aggregated) is difficult and is not only depending on their category, but also (and maybe this is even more important) on their location (Kuhn, 2012). (As a matter of fact, atoms will not be aggregated into objects, because all links, including aggregation, only exists between objects. Neither in the coverage-based view nor in the object-based view atoms "live" on their own. They are either embedded into a grid or into an object. This means, initially, every atom makes up an object). Linking depends on the application/domain and is defined by a concept of **neighbourhood**. This might be simply a 4-connected neighbourhood of a grid or a more complicated concept.

## How to use it?

In the best case you can consider it as "experimental". It is currently not operational and the main intention (for now) is also not to make it operational. The main reason is that it is very slow since it more or less directly implements the GIScience theory without considering performance aspects.

### Usage

If you want to use it, the following environment is required:

Interpreter: 

```bash
python >3.5
```

Packages:

```bash
gdal
numpy
shapely
scikit-image
```

Clone the repository by using the following command:

```bash
git clone git@github.com:whisperingpixel/rhino.git
```

And set up your environment on conda or whatever you want to have. Optionally, you can run the tests:

```bash
python ./rhino-test.py
```

Try out the demo, which is provided in the repository. See the result, which is generated:

```
python ./demo.py
```



### Structure

RHINO is completely object-oriented. The main classes are in the file *rhino.py*, the files *rhino_tools.py* contain additional procedures, which are not directly related to what RHINO does and are helpers (such as printing a nice progress bar etc).

Everything of the aforementioned type Datacube, Atom, Object, Coverage, Link, Neighbourhood is a class in rhino.py

## Is it possible to share, have or steal it?

This is the core of my PhD ... input is welcome, but please don't steal it or at least wait until I got my PhD :-)

## References

Blaschke, T., 2010. Object based image analysis for remote sensing. *ISPRS journal of photogrammetry and remote sensing*, *65*(1), pp.2-16. https://www.sciencedirect.com/science/article/pii/S0924271609000884

Goodchild, M. F., Yuan, M., & Cova, T. J. (2007). Towards a general theory of geographic representation in GIS. *International journal of geographical information science*, *21*(3), 239-260. https://www.tandfonline.com/doi/abs/10.1080/13658810600965271 

Kuhn, W., 2012. Core concepts of spatial information for transdisciplinary research. *International Journal of Geographical Information Science*, *26*(12), pp.2267-2276. https://www.tandfonline.com/doi/full/10.1080/13658816.2012.722637

Lang, S., Kienberger, S., Tiede, D., Hagenlocher, M. and Pernkopf, L., 2014. Geons–domain-specific regionalization of space. *Cartography and Geographic Information Science*, *41*(3), pp.214-226. https://www.tandfonline.com/doi/abs/10.1080/15230406.2014.902755

# Contact

Please write me if you have any questions or encounter any problems with RHINO.

```
Mr. Martin Sudmanns
University of Salzburg
Department of Geoinformatics - Z_GIS
Integrated Spatial Analysis
martin.sudmanns@sbg.ac.at
```

