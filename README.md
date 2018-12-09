[Click here for the full description][home]

This file is part of the geometric calculation library.

The geometric calculation library uses David Avis' [LRS library][lrs], 
to compute the vertex enumeration of a polyhedral set P = {x: A*x<=b}. 
It also computes the facet enumeration to P = conv(V) + cone(R). 
It has a function to compute the projection of a polyhedron, as well as reduction functions
to obtain minimal representations of both V and H-representations of
polyhedra.

INSTALL:
The geometric calculation library offers an interface to Matlab for all
its functionality. It requires the GMP library.

Create a text file named User.make containing the variables MATLABROOT and INSTALLDIR
(optionally you might need GMPINCLUDE) which provide the absolute path to your Matlab 
installation and to the desired installation directory respectively (and the location of 
the GMP header `gmp.h`).

Run 
```	
	make
```

You should now be able to run the mex files from Matlab.

[lrs]: http://cgm.cs.mcgill.ca/~avis/C/lrs.html
[home]: http://worc4021.github.io/GeoCalcLib