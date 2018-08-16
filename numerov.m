%-------------------------------------------------------------------------
% Function name : numerov
% Description   : giving the result of second order differential equation
% which the first order term does not appear.
% Input         : y0: first boundary condition in differential equation in
% which it is the initial condition. y1 : first iteration of the y, this is
% sometimes set as y0
% output        : psi = the array of solution of numerov algorithm
% author        : Norma Dani Risdiandita
% date          : 20/03/2018
%
% Numerov Algorithm in general could be defined as the formula
%
%  y(x+h) = 2*(1-(5/12)*h^2*k^2(x))*y(x) - (1+(h^2/12)*k^2*(x-h))*y(x-h)
%           ____________________________________________________________
% 
%                                 1+h^2*k^2(x+h)/12
%
% as a solution of second order ordinary differential equation
%                d^2y = k^2 y
%                ____   
%                dx^2
%-------------------------------------------------------------------------

function [psi] = numerov(y0,y1,x,f,h)
    psiinitial = zeros(1,length(x));
    psiinitial(1) = y0;
    psiinitial(2) = y1;
    for n = 2:length(x)-1
        left = 2.*(1-(5/12).*h*h*f(x(n)));
        right = 1+(1/12).*h*h*f(x(n-1));
        down = 1 + (1/12).*h.*h.*f(n+1);
        psiinitial(n+1) = (left.*psiinitial(n)-right.*psiinitial(n-1))/down;
    end
    psi = psiinitial;
end
