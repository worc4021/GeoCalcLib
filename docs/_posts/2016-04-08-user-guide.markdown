---
layout: post
title:  "User Guide"
date: 2016-04-09 12:32:35 +0100
categories: jekyll update
---

> Warning: There are no guarantees for the robustness of the GeoCalcLib!
> Please report any malfunction.

The GeoCalcLib provides five different wrappers to call the LRS routines:

1. `facetEnumeration()` calculates a half space description for the provides vertex/ray description.

2. `vertexEnumeration()`calculates a vertex description for the provided half space description.

3. `inequalityReduction()` produces an irredundant half space description for a polyhedron in half space description.
4. `vertexReduction()` produces an irredundant vertex/ray description of a polyhedron in vertex/ray description.
5. `projectPolyhedron()` does not call a LRS routine, but projects a polyhedron in half space description and returns the result in half space description.
6. `LRS()` allows passing more data to the LRS engine to perform the vertex and facet enumeration. Allows passing parameters.

On top of the direct interface to the LRS library `createLRSfile()` and `readLRSfile()` are provided to write and read an LRS input file respectively.


Facet Enumeration
==================

Generate a set of vertices:

\begin{equation}
V = \left(\begin{array}{c}
	v_1^T \\\\[0em]
	v_2^T \\\\[0em]
	\vdots
	\end{array}\right)
\end{equation}

For example $$P = \text{conv}\left\{\left(\begin{array}{c} 1 \\\\[0em] 0 \end{array}\right),\left(\begin{array}{c} 0 \\\\[0em] 1 \end{array}\right), \left(\begin{array}{c} -1 \\\\[0em] 0 \end{array}\right), \left(\begin{array}{c} 0 \\\\[0em] -1 \end{array}\right)\right\}$$ translates to

{% highlight Matlab %}
    V = [eye(2);-eye(2)];
{% endhighlight %}

We can produce a plot of the polytope $$P$$ using 

{% highlight Matlab %}
	k = convhull(V(:,1),V(:,2));
	fill(V(k,1),V(k,2),'r');
{% endhighlight %}


{% include raute.svg %}

The H-representation of $$P$$ is trivially obtained as

\begin{equation}
P = \left\\{x\in\mathbb R^2 : \begin{pmatrix}
1 & 1 \\\\[0em]
-1 & 1 \\\\[0em]
1 & -1 \\\\[0em]
-1 & -1
\end{pmatrix} x \leq \begin{pmatrix} 1 \\\\[0em] 1 \\\\[0em] 1 \\\\[0em] 1 \end{pmatrix} \right\\}
\end{equation}

This can result is obtained by calling `[A,b] = facetEnumeration(V)`. Called with one argument `facetEnumeration(V)` assumes that all rows in `V` are vertices. 

If rays are present, an additional vector `type` must be provided 

{% highlight Matlab %}
[A,b] = facetEnumeration(V,type)
{% endhighlight %}


each row in `type` is either `type(i) = 1` if `V(i,:)` is a vertex or `type(i) = 0` if `V(i,:)` is a ray.
`[A,b] = facetEnumeration(V,ones(size(V,1),1))` produces the same result as `[A,b]=facetEnumeration(V)`.


Vertex Enumeration
==================

The `vertexEnumeration()` function produces a vertex/ray representation for a given polyhedron $$P = \left\{x:Ax\leq b\right\}$$.

Assume we want to understand the epigraph $$\text{epi}(\varphi)$$ of the piecewise affine function $$\varphi(x) = \max\{x+1,2x,3x-4\}$$
for $$x\in[0,5]$$.
The epigraph is given by

\begin{equation}
\text{epi}(\varphi) = \left\\{ (x,t) : 0\leq x\leq 5, x+1\leq t, 2x \leq t, 3x-4\leq t
\right\\}
\end{equation}

\begin{equation}
A = \begin{pmatrix}
	1 &   -1\\\\[0em]
     2 &   -1\\\\[0em]
     3 &   -1\\\\[0em]
     1  &   0\\\\[0em]
    -1  &   0
\end{pmatrix}
, b = \begin{pmatrix}-1 \\\\[0em] 0 \\\\[0em] 4 \\\\[0em] 5 \\\\[0em] 0
\end{pmatrix}
\end{equation}

Using this as an input

{% highlight Matlab %}
>>[V,type] = vertexEnumeration(A,b)
V =

     0     1
     0     1
     1     2
     4     8
     5    11
     0     1

type =

     1
     0
     1
     1
     1
     0
{% endhighlight %}

we obtain the vertices and rays generating the epigraph. 

{% highlight Matlab %}
plot(V(logical(t),1),V(logical(t),2))
{% endhighlight %}

{% include epigraph.svg %}

Inequality Reduction
====================

The call `[Aout,bout] = inequalityReduction(Ain,bin)` returns an irredundant H-representation of $$P = \{x: A_{in}x\leq b_{in}\}$$.

Vertex Reduction
================

The call `[Vout,tout] = vertexReduction(Vin,tin)` returns an irredundant V-representation of $$ P = \text{conv}\{v_{i}\}\oplus \text{cone}\{r_{i}\}$$ where the vertices in `Vin` are passed by setting `tin(i) = 1` if `Vin(i,:)` is a vertex and `tin(i) = 0` if `Vin(i,:)` is a ray.

If $$ P = \text{conv}\{v_i\} $$, passing `tin = ones(size(V,1),1)` can be omitted `Vout = vertexReduction(Vin)` assumes all passed points are vertices.


Projection
==========

The projection of polyhedral sets is performed for polyhedra in H-representation $$P = \{x\in\mathbb R^n: A_{in}x \leq b_{in}\}$$,
and is always projected on to $$\mathbb R^{n-d}$$. The projected set $$\pi_d(P) = \{y\in\mathbb R^{n-d}: \exists z\in\mathbb R^d y\times z\in P\}$$ is returned in H-representation.

Internally a vertex enumeration is performed on $$P$$, the projection operator $$\pi_d$$ acts trivially on points by just dropping off the last $$n-d$$ elements, after that a facet enumeration is appended to produce the H-representation of $$\pi_d(P)$$.
The function call `[Aout,bout] = projectPolyhedron(Ain,bin,d)` performs this operation for general polyhedra:

{% highlight Matlab %}
P = gallery('uniformdata',[30,3],1);
[A,b] = facetEnumeration(P);

K = convhull(P);
figure(1)
trisurf(K,P(:,1),P(:,2),P(:,3))
hold on

[C,d] = projectPolyhedron(A,b,1);

V = vertexEnumeration(C,d);
k = convhull(V(:,1),V(:,2));
fill(V(k,1),V(k,2),'b');
hold off
{% endhighlight %}


{% include projection.svg %}


Vertex and facet enumeration with LRS parameters
================================================

The functions `vertexEnumeration()` and `facetEnumeration()` perform a vertex and facet enumeration respectively for data that was preconditioned by the user.
The `vertexEnumeration()` call requires the user to prepare a type vector if rays are present etc. and the `facetEnumeration()` call only treats inequality constrained sets.
To exploit (almost) the full scope of possibilities the LRS engine offers the `LRS()` function is provided.

The argument is a `struct` with the required field `'rep'` specifying whether a V- or an H-representation is passed.
Each representation has its own arguments:


Assume we want to facet enumerate the first quadrant, i.e. $$P=\text{cone}\left\{\begin{pmatrix}1 \\\\ 0\end{pmatrix}, \begin{pmatrix}0 \\\\ 1\end{pmatrix} \right\}$$, then we define the structure

{% highlight matlab %}
s = struct('rep','V','R',eye(2));
[A,b] = LRS(s)
A =

     0    -1
    -1     0


b =

     0
     0
{% endhighlight %}

If we have $$P = \text{cone}\left\{\begin{pmatrix}1 \\\\ 0\end{pmatrix}, \begin{pmatrix}0 \\\\ 1\end{pmatrix} \right\} \oplus \text{conv}\left\{\begin{pmatrix}1\\\\1 \end{pmatrix} \right\}$$ we set the additional field 

{% highlight matlab %}
s.V = [1,1];
[A,b] = LRS(s)
A =

     0     0
    -1     0
     0    -1


b =

     1
    -1
    -1
{% endhighlight %}

That is, if the `s.rep='V'` at least one of `s.V` or `s.R` has to be passed.


Assume now that we want to vertex enumerate the 3-dimensional simplex $$P = \{x:x_1\geq0,x_2\geq0,x_3\geq0,x_1+x_2+x_3=1\}$$,
that is the representation of `rep='H'`

{% highlight matlab %}
s = struct('rep','H','Aineq',-eye(3),'bineq',zeros(3,1),'Aeq',ones(1,3),'beq',1);
[V,t] = LRS(s)
V =

     1     0     0
     0     1     0
     0     0     1


t =

     1
     1
     1
{% endhighlight %}

The interpretation of the output for a vertex enumeration is exactly the same as for `vertexEnumeration()` that is, `V(i,:)` is a vertex if `t(i)=1` and `V(i,:)` is a ray if `t(i)=0`.

If `s.rep='H'` then `s.Aineq` and `s.bineq` are required, `s.Aeq` and `s.beq` are optional.


In addition to the computational data parameters may be passed, so far `maxcobases`, `maxoutput` and `maxdepth` are accepted. See [LRS options][lrs-options] for details.

If we want to restrict the number of search depth in the vertex enumeration of a 12-dimensional cube, and further we also want to constrain the number of explored cobases and only need 5 vertices we would call

{% highlight matlab %}
[V,t] = LRS(struct('rep','H','Aineq',[eye(12);-eye(12)],'bineq',ones(24,1),'maxdepth',3,'maxcobases',12,'maxoutput',5))

V =

    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1
     1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1
    -1     1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1
     1     1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1
    -1    -1     1    -1    -1    -1    -1    -1    -1    -1    -1    -1


t =

     1
     1
     1
     1
     1
{% endhighlight %}


> Warning: Although this function has been thoroughly tested, there might be uncaught exceptions. 
> Matlab will immediately crash if you come across an uncaught exception.
> These exceptions have no numerical nature, but can occur when passing unaccepted data types *(passing cell arrays instead of arrays etc.*).
> Please test your code your function calls before calling this function in an automated environment!


Writing and reading LRS input files
===================================

Many LRS users might prefer to call the LRS executable directly in order to be able to see the output as it is obtained, rather than waiting for the entire solution to be returned. 
For this purpose `createLRSfile('param1',param1,...)` is provided to generate a LRS input file from within the Matlab environment.

The supported parameters are

- `fname` - specifies the file name of the LRS input file created, if omitted a file called `lrstest.ine` is created in the current directory.
- `rep` - `'V'` or `'H'` implying that the data is passed in a V- or H-representation respectively.
- If `rep='V'` vertices are passed as the parameter `V` in a matrix `vdat` containing a vertex/ray in each row, optionally `Type` can be used to pass a vector of `1/0` such that `type(i) = 0` implies `vdat(i,:)` is a ray and `type(i) = 1` implies `vdat(i,:)` is a vertex.
If `Type`is not specified all rows of `vdat` are assumed to be vertices.
- If `rep='H'` the set of inequalities $$\{x:A_{dat}x\leq b_{dat}\}$$ is passed by setting `A` to `Adat` and `b` to `bdat`.
If `b` is omitted the set is assumed to be $$\{x:A_{dat}x\leq{\bf{1}}\}$$
- The optional parameter `tol` may be set to the desired conversion tolerance to a rational format, default `1e-12`.


To create an LRS input file `myfile.ine` to vertex enumerate the $$\ell_\infty$$ unit ball in two dimensions $$\|x\|_\infty\leq 1$$ we would use

{% highlight matlab %}
A = [eye(2);-eye(2)];
b = ones(4,1);
createLRSfile('rep','H','A',A,'b',b,'fname','myfile.ine');
{% endhighlight %}

Users who have produced results using the LRS executable wanting to import it into the Matlab environment would first change the head of the result file into an admissible LRS input file, i.e replace the `***** n rational` line by `m n rational`, where `m` and `n` are the number of rows and columns respectively. Then pass the `filename` to `readLRSfile('filename')`.


To read in the result of a vertex enumeration in `myfile.ext`, i.e. the result is in V-representation, we would call

{% highlight matlab %}
[V,type] = readLRSfile('myfile.ext');
{% endhighlight %}

[lrs-options]: http://cgm.cs.mcgill.ca/%7Eavis/C/lrslib/USERGUIDE.html#Options