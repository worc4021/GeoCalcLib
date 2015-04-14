function [V,type] = vertexEnumeration(A,b)


temp = calllib('libgeocalc','vertexEnumeration',A,b);
V = temp{1};
type = temp{2};