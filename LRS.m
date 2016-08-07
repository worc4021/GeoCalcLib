function [Data,data,volume] = LRS(s)
% [Data,data,volume] = LRS(s)
% Direct LRS call.
% Structure s is used to pass data to the LRS engine.
%
% For V-representation:
% s.rep = 'V'
% s.V = V Matrix holding all vertices, each row is one vertex.
% s.R = R Matrix holding all rays, each row is one ray.
% At least one of the two must be non-trivial.
% Result is H-representation {x:Data*x<=data}.
% If s.getvolume is passed volume holds the volume of the polytope.
%
% For H-representation:
% s.rep = 'H'
% s.Aineq and s.bineq represent {x:s.Aineq*x<=s.bineq}.
% s.Aeq and s.beq represent {x:s.Aeq*x==beq}.
% The inequalities must not be omitted, the equalities are optional.
% Result is V-representation Data(i,:) is vertex if data(i)==1, ray
% for data(i) == 0.
% 
% Additional parameters are:
% s.maxdepth, specifies trancation of treesearch.
% s.maxoutput, truncates maximal number of output.
% s.maxcobases, limits the number of explored cobases.


temp = calllib('libgeocalc','lrsCall',s);
Data = temp{1};
data = temp{2};
volume = temp{3};