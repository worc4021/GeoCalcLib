function [Vout,typeOut] = vertexReduction(varargin)
% [Vout,typeOut] = vertexReduction(Vin,typeIn)
% [Vout,typeOut] = vertexReduction(Vin)
% If no typeIn is passed all rows of Vin are assumed to be vertices.
% Otherwise each row in Vin(i,:) corresponds to a vertex if typeIn(i) = 1,
% or a ray if typeIn(i) = 0. In Vout every row corresponds to either a 
% vertex or a row, where in typeOut 1 denotes the type vertex and 0 
% denotes the type ray.

V = varargin{1};
if nargin <2
    type = ones(size(V,1),1);
else
    type = varargin{2};
end

temp = calllib('libgeocalc','vertexReduction',V,type);

Vout = temp{1};
typeOut = temp{2};