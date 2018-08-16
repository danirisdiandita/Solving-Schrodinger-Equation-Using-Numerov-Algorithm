%-------------------------------------------------------------------------
% Function name : trianglepotential
% Description   : a function of triangle wave
% Input         : x = an array; a = period of the triangle wave
% output        : output = a triangle function as an array
% author        : Norma Dani Risdiandita
% date          : 20/03/2018
% the triangle potential is in the unit of hbar^2/m
%-------------------------------------------------------------------------
function output = trianglepotential(a,x)
    output = (2/a).*(x./a - floor(x./a + 0.5)).*(-1).^floor(x./a + 0.5);
end

