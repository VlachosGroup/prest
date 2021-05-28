Reaction Mechanism Generator - Python Reaction Stencil (pReSt)
==============================================================

The **p**ython **Re**action **St**encil (pReSt) generates the reaction rules using a set of elementary steps for a reaction system. 

Documentation
--------------

See our [documentation page](https://github.com/VlachosGroup/pReSt/wiki/python-reaction-stencil-(pReSt)-Usage-Instructions) for examples.

Developer
---------
Udit Gupta (ugupta@udel.edu)

Dependencies
------------

- Python3
- [ASE](https://wiki.fysik.dtu.dk/ase/about.html) : Used for reading in 3-D molecular structures
- [Numpy](http://www.numpy.org/) : Used for vector and matrix operations
- [rdkit](http://www.rdkit.org/) : Used for performing bond transformations
- [NetworkX](https://networkx.org/) : Used for performing graph manipulation

Getting Started
---------------
1. Install using pip::

  >  pip install pReSt
 
Files for using the Visualization tool
--------------------------------------
1) species.xyz - Species files containing 3-D coordinates for each molecule
2) reaction.txt - A list of elementary steps containing species provided along with species.xyz file

Features in the Reaction Mechanism Generator:
---------------------------------------------
- Reaction rule generation
- Reaction network generation using reaction rules from a database
- Checking if a reaction exists in the database
- Generating 3-D molecular structures using SMILES strings (**In Progress**)

License
-------

This project is licensed under the GNU LGPL License - see the [LICENSE.md](https://github.com/VlachosGroup/pReSt/blob/master/LICENSE.md) file for details

Contributing
------------

If you have a suggestion or find a bug, please post to our Issues page with 
the enhancement or bug tag respectively.

Finally, if you would like to add to the body of code, please:

- fork the development branch
- make the desired changes
- write the appropriate unit tests
- submit a pull request.


Questions
---------

If you are having issues, please post to our Issues page with the 
help wanted or question tag. We will do our best to assist.

Funding
-------

This material is based upon work supported by the Department of Energy's Office 
of Energy Efficient and Renewable Energy's Advanced Manufacturing Office under 
Award Number DE-EE0007888-9.5.

Special Thanks
--------------

-  Guen Ho Gu (Discussion)
-  Qiang Li (Discussion)

