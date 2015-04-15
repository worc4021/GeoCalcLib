function [Aout,bout] = projectPolyhedron(A,b,dim)
% [Aout,bout] = projectPolyhedron(A,b,dim)
% Returns the projection of the n dimensional Polyhedron P = {x: A*x<=b}  
% onto the n-dim dimensional Polyhedron pP = {x: Aout*x<=bout}.

temp = calllib('libgeocalc','projection',A,b,uint32(dim));

Aout = temp{1};
bout = temp{2};