function [A,b] = facetEnumeration(varargin)

V = varargin{1};

if nargin<2
    type = ones(size(varargin{1},1),1);
else
    type = varargin{2};
end

temp = calllib('libgeocalc','facetEnumeration',V,type);
    
A = temp{1};
b = temp{2};