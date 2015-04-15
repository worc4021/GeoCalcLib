function [A,b] = facetEnumeration(varargin)
% [Aout,bout] = facetEnumeration(V)
% [Aout,bout] = facetEnumeration(V,type)
% Compute the facet enumeration of conv_i(V(i,:)') where all rows
% in V are taken as vertices. If 'type' is passed it has to describe
% the type of the associated row passed in V, 1 corresponds to a vertex,
% and 0 corresponds to a ray.


V = varargin{1};

if nargin<2
    type = ones(size(varargin{1},1),1);
else
    type = varargin{2};
end

temp = calllib('libgeocalc','facetEnumeration',V,type);
    
A = temp{1};
b = temp{2};