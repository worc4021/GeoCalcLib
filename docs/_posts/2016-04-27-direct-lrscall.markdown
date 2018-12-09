---
layout: post
title:  "Example"
date: 2016-04-27 17:31:35 +0100
categories: jekyll update
---

> Warning: There are no guarantees for the robustness of the GeoCalcLib!
> Please report any malfunction.

The permutahedron[^1] $$\Pi_{d-1}\subseteq\mathbb R^d$$ is given by the convex hull of all permutations of the coordinates of 

\begin{equation}
\begin{pmatrix} 1 \\\\ 2 \\\\ \vdots \\\\ d \end{pmatrix}
\end{equation}

It is easy to see that all such vectors satisfy 

\begin{equation}
\begin{pmatrix} 1 && 1 && \cdots && 1 \end{pmatrix}x = \frac{d(d-1)}{2}
\end{equation}

i.e. $$\Pi_{d-1}\subseteq \{x\in\mathbb R^d: {\bf{1}}x= \frac{d(d-1)}{2}\}$$ and therefore has a dimension of $$d-1$$.

We use $$\Pi_3$$ as an example. Naturally $$\Pi_3$$ is a 3-dimensional polytope embedded in $$\mathbb R^4$$, in order to visualise it we need to project it onto $$\{x\in\mathbb R^d: {\bf{1}}x= \frac{d(d-1)}{2}\}$$. This can easily be done using a QR decomposition of $${\bf{1}}^T$$. If $$QR={\bf{1}}^T$$ then $$Q^Tv_i$$ rotates $$v_i$$ in such a way that its first coordinate corresponds to $$c\frac{d(d-1)}{2}$$. We can drop the first coordinate and study at the convex hull of 3-dimensional vectors.

In Matlab we can use

{% highlight matlab %}
V = [kron((1:4)',ones(64,1)),repmat(kron((1:4)',ones(16,1)),4,1),repmat(kron((1:4)',ones(4,1)),16,1),repmat(kron(ones(4,1),(1:4)'),16,1)];

idx = false(256,1);
for i = 1:256
    if length(unique(V(i,:))) == 4
        idx(i) = true;
    end
end

V = V(idx,:);

[Q,~] = qr([1,1,1,1]');

V = V*Q;
V = V(:,2:4);

K = convhull(V);
trisurf(K,V(:,1),V(:,2),V(:,3))
{% endhighlight %}

{% include permutahedron.svg %}

Assume now we want to perform a projective transformation, that is we want to intersect the _homogenisation_ of $$\Pi_3$$ with a hyperplane $$H = \{x\in\mathbb R^4 : ax = b\}$$. Since $$\text{homog}(\Pi_3)$$ is a cone in $$\mathbb R^4$$ simply rotating $$\text{vert}(\text{homog}(\Pi_3))$$ using a QR decomposition of $$a^T$$ and discarding the first coordinates would yield $$\mathbb R^3$$.

Therefore we use a facet enumeration of $$\Pi_3$$ to produce $$\text{homog}(\Pi_3)$$, with that we will enumerate the vertices of 

\begin{equation}
T(\Pi_3,H) = \\{x\in\mathbb R^4: x\in\text{homog}(\Pi_3)\wedge x\in H\\}
\end{equation}

Naturally $$\text{vert}(T(\Pi_3,H))\in H$$ such that we can again project them by using a QR decomposition.

For the numerical example assume $$H = \{\begin{pmatrix} 1 & 2 & 3 &15 \end{pmatrix}x=1\}$$. We use the following Matlab code:

{% highlight matlab %}
s.rep = 'V';
s.V = V;
[PiA,Pib] = LRS(s);

homogA = [zeros(1,3),-1;...
          PiA,-Pib];
homogb = zeros(size(homogA,1),1);

clear s

s.rep = 'H';
s.Aineq = homogA;
s.bineq = homogb;
s.Aeq = [1,2,3,15];
s.beq = 1;
[TV,t] = LRS(s);

[Q,~] = qr(s.Aeq');
TV = TV*Q;
TV = TV(:,2:4);

K = convhull(TV);
trisurf(K,TV(:,1),TV(:,2),TV(:,3))
{% endhighlight %}

{% include permutahedronTrans.svg %}


[^1]: G. M. Ziegler - _Lectures on Polytopes_