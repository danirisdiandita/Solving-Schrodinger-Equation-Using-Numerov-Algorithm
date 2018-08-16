%-------------------------------------------------------------------------
% Function name : normalize
% Description   : giving the normalized array from the input of an array
% Input         : x = an array
% output        : output = a normalized array.
% author        : Norma Dani Risdiandita
% date          : 20/03/2018
%-------------------------------------------------------------------------
function [output] = normalize(x)
    constant = dot(x,x);
    output = x/sqrt(constant);
end

