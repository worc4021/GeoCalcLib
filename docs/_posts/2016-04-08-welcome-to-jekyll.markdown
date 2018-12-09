---
layout: post
title:  "Installing the GeoCalcLib"
date:   2016-04-08 12:32:35 +0100
categories: jekyll update
---

<!-- >[!WARNING]
> This way of setting up the GeoCalcLib has been deprecated due to maintanence difficulties.


> Warning: There are no guarantees for the robustness of the GeoCalcLib!
> Please report any malfunction.

The [LRS library][lrs-link] works on both Mac and Linux and can be made to work on Windows.

In order to work you will need the [GMP Library][gmp] installed on your system. 
During the build of the GMP library, make sure you install it with the same ABI as your Matlab, i.e. if you have a 64-bit Matlab make sure you select `ABI=64` during your build. 
If you are unsure which Matlab version you have, run -->

<!-- {% highlight matlab %}
>>ver matlab
{% endhighlight %}

and check which Java version is linked, they are installed compatibly. The GeoCalcLib is untested on 32-bit installations. -->

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
MATLABROOT = /Applications/MATLAB_R2018b.app

# Path to which everything should be installed, has to be on Matlab path!
INSTALLDIR = /Users/Username/Documents/MATLAB/Functions
{% endhighlight %}

If necessary add `GMPINCLUDE` with the path to the `gmp.h` file (I only had to do this on osx).

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

<!-- Now you can run 

{% highlight bash %}
make
{% endhighlight %}

to build the interface.

Next you have to load the library in Matlab, for this you need to specify the location of the header file containing the information Matlab needs to call the individual functions, see [here][matlab-dll] for details.
You can load the interface on startup, to do this run Matlab and edit the startup.m file:

{% highlight matlab %}
>> edit startup
{% endhighlight %}

If it already exists add the lines

{% highlight matlab %}
[~,~] = loadlibrary('libgeocalc','PATH/mainFunctions.h');
display('GeoCalc library loaded.');
{% endhighlight %}

where you replace `PATH` with the absolute path to the directory containing the header file `mainFunctions.h`.

Once the library is loaded you need to unload it before you can exit Matlab! So exit Matlab from the session in which you modified `startup.m` *(it has not been executed, therefore GeoCalcLib has not been loaded*). 

Run Matlab again and edit `finish.m` by typing

{% highlight matlab %}
edit finish
{% endhighlight %}

add the line 

{% highlight matlab %}
unloadlibrary('libgeocalc');
{% endhighlight %}


If you do not wish to load the GeoCalcLib every time you run Matlab, load the library by calling 

{% highlight matlab %}
[~,~] = loadlibrary('libgeocalc','PATH/mainFunctions.h');
{% endhighlight %}

at the beginning of your session and 

{% highlight matlab %}
unloadlibrary('libgeocalc');
{% endhighlight %}

at the end of your session. -->


[lrs-link]: http://cgm.cs.mcgill.ca/~avis/C/lrs.html
[matlab-dll]: http://uk.mathworks.com/help/matlab/using-c-shared-library-functions-in-matlab-.html
[gmp]: https://gmplib.org
[matlab-path]: http://uk.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html