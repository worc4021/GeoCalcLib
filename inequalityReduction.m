function [Aout,bout] = inequalityReduction(A,b)
% [Aout,bout] = inequalityReduction(Ain,bin)
% Returns the minimal representation of {x:Ain*x<=bin}.

temp = calllib('libgeocalc','ineqReduction',A,b);

Aout = temp{1};
bout = temp{2};