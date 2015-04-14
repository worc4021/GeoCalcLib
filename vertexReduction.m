function [Vout,typeOut] = vertexReduction(varargin)

V = varargin{1};
if nargin <2
    type = ones(size(V,1),1);
else
    type = varargin{2};
end

temp = calllib('libgeocalc','vertexReduction',V,type);

Vout = temp{1};
typeOut = temp{2};