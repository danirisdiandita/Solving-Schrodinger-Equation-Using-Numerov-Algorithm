%-------------------------------------------------------------------------
% Function name : sawtoothpotential
% Description   : a function of sawtooth wave
% Input         : x = an array; a = period of the sawtooth wave
% output        : output = a sawtooth function as an array
% author        : Norma Dani Risdiandita
% date          : 20/03/2018
% the potential is in the unit of h^2/m
%-------------------------------------------------------------------------
function output = sawtoothpotential(a,x)

    output = 1.*(x./a - floor(0.5 + x./a));
    
end