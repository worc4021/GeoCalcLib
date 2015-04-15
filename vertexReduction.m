function [Vout,typeOut] = vertexReduction(varargin)
% [Vout,typeOut] = vertexReduction(Ain,bin)
% [Vout,typeOut] = vertexReduction(Ain)
% In no bin is passed the constant one vector is assumed.
% In Vout every row corresponds to either a vertex or a row
% of {x:Ain*x<=bin}, where in typeOut 1 denotes the type
% vertex and 0 denotes the type ray.

V = varargin{1};
if nargin <2
    type = ones(size(V,1),1);
else
    type = varargin{2};
end

temp = calllib('libgeocalc','vertexReduction',V,type);

Vout = temp{1};
typeOut = temp{2};