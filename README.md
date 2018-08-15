# RHINO

## What is RHINO?

Rhino is a translation layer, which allows switching between a coverage view and an object view within a data cube.

## What does RHINO mean?

The name RHINO does not mean anything (or everything, depending on your attitude and whether you are more optimistic or pessimistic).

## Why do I want to have this?

Geographers (and people from other domains, too) traditionally separate their view of the world either into a coverage-based or object-based representation (Sometimes, this is referred to as raster and vector, although it is not quite the same). This is documented by a number of software and approaches, which support usually only one of them. However, there are several application areas, where an analyst has to switch between both worlds seamlessly. RHINO solves exactly this gap by building a bridge between these two worlds: It does not matter whether you want to perform object-based or coverage-based analysis. In RHINO you can conduct both analyses and can even combine them into an integrated workflow. 

Now, for the sake of completeness, it is true that modern GIS or remote-sensing software (ArcGIS, eCognition) allows switching between the two worlds to some degree. However, RHINO is especially designed to be built on top of data cubes and implements some of the more recent considerations about the nature of geospatial information coming from GIScience theory.

## How does it work?

The most basic element in rhino is an **atom**. An atom relates a location on the Earth to a property and its value. For example at location [lat] / [lon] the property *land cover* might be *water*. This is what an atom can tell. Both of the two representations of the world usually consist of a lot of atoms, but it depends on your preferred representation of the world what you do with them.

You can put all of the the atoms into a regular or irregular grid, which you can treat as a synonym of a coverage for now.  You chose to view the atoms in a **coverage**-based representation. 

On the other hand, you can also put all of the surrounding atoms of the aforementioned water atom, which have the same property value, into a common "container" and draw a border around them. Now you chose to have an **object**-based representation. An object has properties such as size, shape, etc.

Having an objects also means to have relationships (associations) between them. In RHINO they are called **link**. The reason is that linking objects might generate a new object (e.g. by generalisation of surface water bodies) or maintain a relationship between them (e.g. a contract between a cloud and cloud shadow).

Determining which objects are aggregated or generalised (as a matter of fact, atoms will not be aggregated into objects, because aggregation only exists between objects. Neither in the coverage-based view nor in the object-based view atoms "live" on their own. They are either embedded into a grid or into an object. This means, initially, every atom makes up an object) is difficult, because it depends on the application/domain and is defined by a **neighbourhood**-concept. This might be simply a 4-connected neighbourhood of a grid or a more complicated concept.

