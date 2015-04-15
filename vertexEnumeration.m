function [V,type] = vertexEnumeration(A,b)
% Returns all vertices/rays of {x:A*x<=b} in V row wise,
% the type decodes vertices with 1 and rays with 0.


temp = calllib('libgeocalc','vertexEnumeration',A,b);
V = temp{1};
type = temp{2};