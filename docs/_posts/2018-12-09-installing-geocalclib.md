---
layout: post
title:  "Installing the GeoCalcLib - Mex approach"
date:   2018-12-08 12:32:35 +0100
categories: install update
---


> Warning: The previous way of installing the GeoCalcLib has been deprecated. This is the maintaned method of building the binaries.


A prerequisite to run the GeoCalcLib is to have the [GMP Library][gmp] on your machine. At build time, you will need the main header file `gmp.h` however the library itself is largely irrelevant, as Matlab provides a gmp redistributable.

To install the GeoCalcLib open your Terminal navigate to your desired directory and clone the repository by typing

{% highlight bash %}
git clone https://github.com/worc4021/GeoCalcLib.git
cd GeoCalcLib
{% endhighlight %}

Create a text file named `User.make` in which you define the absolute path to your Matlab installation `MATLABROOT`
and the desired installation directory `INSTALLDIR` *(a directory on your Matlab path*), e.g.

{% highlight make %}
# Specify the absolute path to your Matlab installation
MATLABROOT = /Applications/MATLAB_R2019b.app

# Path to which everything should be installed, has to be on Matlab path!
INSTALLDIR = /Users/Username/Documents/MATLAB/Functions
{% endhighlight %}

If necessary add `GMPINCLUDE` with the path to the `gmp.h` file (I only had to do this on osx ).

Run make

{% highlight bash %}
make
{% endhighlight %}

If you do not yet have a directory on your Matlab path to store your functions and want to add one, run `edit startup` in Matlab and add 

{%highlight matlab %}
addpath('/absolute/path/to/your/directory');
{% endhighlight %}

or find [alternative ways of setting up your path here][matlab-path].

You should now be set up.


[gmp]: https://gmplib.org
[matlab-path]: http://uk.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html