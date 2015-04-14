function [Aout,bout] = inequalityReduction(A,b)

temp = calllib('libgeocalc','ineqReduction',A,b);

Aout = temp{1};
bout = temp{2};