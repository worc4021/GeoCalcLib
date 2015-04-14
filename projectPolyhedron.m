function [Aout,bout] = projectPolyhedron(A,b,dim)

temp = calllib('libgeocalc','projection',A,b,uint32(dim));

Aout = temp{1};
bout = temp{2};