clear;
% -------------------------------------------------------------------------
% this is a program about solving the Schrodinger Equation of "unusual"
% potential wells using Numerov Algorithm, Shooting Method, and Newton Method.
% The potentials I choose are
% 1. Infinite Square Well
% 2. Infinite Square Well with sawtooth potential
% 3. Infinite Square Well with triangle potential
% The original time independent Schrodinger equation (TISE) is about solving a
% second order ordinary differential equation. The equation is given by
%
% -(hbar^2/2m)*d^2 psi/dx^2 + V(x)*psi = E*psi
% ==> d^2psi/dx^2 = -2*(E-V(x))/hbar^2 = -k^2 *psi ...(1)
%
% which
% hbar = Planck Constant, We set it as 1
% m = mass of a particle, we set it as 1
% psi = the eigenstate of the particle. This parameter is the eigenstate of
% the Hamiltonian of "unusual" potential wells. The solution of 2nd order ODE
% is defined as this. The left and right boundaries of the well has been
% blocked by infinite potential wall. Therefore psi would be zero at left
% and right boundaries as the boundary conditions.
% V(x) = the choosen potential between two Infinite wall. outside the wall
% this value is known to be Infinite. Therefore, there is no possibility that
% the particle will go outside the well.
% E = the eigenvalue of the Hamiltonian. This is the parameter of the Shooting
% for the shooting method.
% -------------------------------------------------------------------------

% DEFINE THE INITIAL PARAMETER
% the width of the well a as the left boundary and b as the right boundary

a = -5.;
b = 5.;

% separating the well coordinate into 1000 points and define the stepsize h,
% and one dimension array x as the array of x coordinate axis.

n = 1000;
h = (b-a)/n;
x = linspace(a,b,n+1);

% Doing the task of Infinite Square Well
%--------------------------------------------------------------------------
% 1. Infinite Square Well
%--------------------------------------------------------------------------
fprintf('1. Infinite Square Well Case\n');

% INITIALIZING THE ARRAY OF OUR RESULT. Numerov Algorithm need to take the
% first two points to make the overall iterations.

psi = zeros(1,n+1);

% INITIALIZING TWO POINTS. psi(1) is the boundary value. Initial psi is 0.
% however the second point is arbitrary since the result from Numerov
% algorithm, an array, need to be normalized to make sure it will
% provide the correct scaling

psi(1) = 0;
psi(2) = 1;

% initializing the energy eigenvalue to be iterated in Shooting Method.
% According to the theory, Energy = Potential + Kinetic.
% since in the Infinite square well, the potential energy is zero between the
% well and the minimum kinetic energy is zero. We can set it as a little bit
% less than zero. Let us say that E = -1

E = -1.;

% the squared k as shown in the equation 1 for Infinite Square Well is
% taken as function handle

k2 = @(x) 2.*E;

% calculating the eigenstate using choosen energy E and normalize the array

psi = numerov(psi(1),psi(2),x,k2,h);
psi = normalize(psi);

% determining the psi at the endpoint as zero.

psiend = 0.0;

% defining the stepsize of the energy for iterating it as E + dE continuously
% until achieving desired result through Shooting method.

dE = 0.01;

% defining the tolerance of psi if it achieves the difference with psiend at
% least less than tolerance. it will return as true result of psi.

tol = 1e-4;

% the 5 eigenstate will be stored in this 5x(length of the psi)

wavefunction1 = zeros(5,length(psi));

% Initializing the array of the 5 lowest energy level as array of zeros.
% the energy levels would be stored in this array.

Earray = zeros(1:5);

% for 5 times Shooting Method

for n = 1:5

    % starting of the Shooting Method.

    while 1
        E = E + dE;
        k2 = @(x) 2.*E;
        psi1 = numerov(psi(1),psi(2),x,k2,h);
        psi1 = normalize(psi1);

        % if the end of psi changed its sign from positive(negative) to
        % negative (positive), the Newton Method will do the stuff by finding
        % the zeros point of psi's endpoint.

        if psi(length(psi)).*psi1(length(psi1)) < 0
            E1 = E - dE;
            E2 = E;
            E1iter = E1;
            E2iter = E2;

            % run the Newton Method until the zero point is found.
            % until the difference of psi's endpoint and psiend less than
            % tolerance factor.

            while abs(psi1(length(psi1))-psiend) > tol
                k2a = @(x) 2.*E1iter;
                k2b = @(x) 2.*E2iter;
                psia = numerov(0,1,x,k2a,h);
                psia = normalize(psia);
                psib = numerov(0,1,x,k2b,h);
                psib = normalize(psib);
                E3iter = E2iter - psib(length(psib)).*(E2iter - E1iter)/(psib(length(psib))-psia(length(psia)));
                k2c = @(x) 2.*E3iter;
                psi1 = numerov(0,1,x,k2c,h);
                psi1 = normalize(psi1);
                E2iter = E3iter;
                E1iter = E2iter;
            end

            %storing the result to psi1. iteration is done by interchanging
            % psi and psi1

            psiresult = psi1;
            fprintf('The energy level at n = %d is %f hbar^2/m \n',n,E3iter);
            Earray(n) = E3iter;
            E = E2;
        end
        k2 = @(x) 2.*E;
        psi = numerov(0.0,1,x,k2,h);
        psi = normalize(psi);
        if abs(psi1(length(psi1))-psiend) < tol
            break;
        end
    end
    %storing to the 5 wavefunction array.
    wavefunction1(n,:) = psiresult;

end

% PLOTTING

figure;
subplot(5,2,1:2:10);
hold on;
plot(x,0.*ones(1,length(x)),'LineWidth',4);
title('First 5 Energy Levels and the Chosen Potential Well');
xlabel('x');
ylabel('Energy (hbar^2/m)');
for m=1:5
    plot(x,Earray(m).*ones(1,length(x)),'LineWidth',4);
end
legend('V_0','E_1','E_2','E_3','E_4','E_5');
hold off;
for n = 1:5
    m = 2.*n;
    subplot(5,2,m);

    % DO NORMALIZATION ACCORDING TO THE QUANTUM NORMALIZATION PROCEDURE AFTER
    % DOING SPLINE INTERPOLATION

    wavefunctionc = @(t) spline(x,wavefunction1(n,:).*wavefunction1(n,:),t);
    wavefunction1(n,:) = wavefunction1(n,:)/sqrt(integral(wavefunctionc,x(1),x(length(x))));
    plot(x,wavefunction1(n,:));
    xlabel('x');
    ytext = strcat('$\psi_',num2str(n),'$');
    ylabel(ytext, 'Interpreter','Latex');
    titletext = strcat('the eigenstate of Infinite Square Well at level ',num2str(n));
    title(titletext);
end

%--------------------------------------------------------------------------
% 2. Infinite Square Well with sawtooth potential
%--------------------------------------------------------------------------
% doing the same method in Infinite Square Well case but with adding sawtooth
% potential to the k2 function handle
fprintf('2. Infinite Square Well Case with sawtooth potential \n');
psi = zeros(1,n+1);
psi(1) = 0;
psi(2) = 1;
E = -1.;

% adding potential to the function handle

k2 = @(x) 2.*(E-sawtoothpotential(2.5,x));
psi = numerov(psi(1),psi(2),x,k2,h);
psi = normalize(psi);
psiend = 0.0;
dE = 0.01;
tol = 1e-4;
wavefunction1 = zeros(5,length(psi));
Earray = zeros(1:10);
for n = 1:5
    while 1
        E = E + dE;
        k2 = @(x) 2.*(E-sawtoothpotential(2.5,x));
        psi1 = numerov(psi(1),psi(2),x,k2,h);
        psi1 = normalize(psi1);
        if psi(length(psi)).*psi1(length(psi1)) < 0
%             E = E - dE;
%             dE = dE/2.;
            % start from here
            E1 = E - dE;
            E2 = E;
            E1iter = E1;
            E2iter = E2;
            while abs(psi1(length(psi1))-psiend) > tol
                k2a = @(x) 2.*(E1iter - sawtoothpotential(2.5,x));
                k2b = @(x) 2.*(E2iter - sawtoothpotential(2.5,x));
                psia = numerov(0,1,x,k2a,h);
                psia = normalize(psia);
                psib = numerov(0,1,x,k2b,h);
                psib = normalize(psib);
                E3iter = E2iter - psib(length(psib)).*(E2iter - E1iter)/(psib(length(psib))-psia(length(psia)));
                k2c = @(x) 2.*(E3iter- sawtoothpotential(2.5,x));
                psi1 = numerov(0,1,x,k2c,h);
                psi1 = normalize(psi1);
                E2iter = E3iter;
                E1iter = E2iter;
            end
            psiresult = psi1;
            fprintf('The energy level at n = %d is %f hbar^2/m \n',n,E3iter);
            Earray(n) = E3iter;
            E = E2;
            % ends here
        end
        k2 = @(x) 2.*(E - sawtoothpotential(2.5,x));
        psi = numerov(0.0,1,x,k2,h);
        psi = normalize(psi);
        if abs(psi1(length(psi1))-psiend) < tol
            break;
        end
    end
    wavefunction1(n,:) = psiresult;
    dE = 0.01;
end

figure;
subplot(5,2,1:2:10);
hold on;
plot(x,sawtoothpotential(2.5,x),'LineWidth',4);
title('First 5 Energy Levels and the Infinite Square Well with sawtooth potential');
xlabel('x');
ylabel('Energy (hbar^2/m)');
for m=1:5
    plot(x,Earray(m).*ones(1,length(x)),'LineWidth',4);
end
legend('V_0','E_1','E_2','E_3','E_4','E_5');
hold off;
for n = 1:5
    m = 2.*n;
    subplot(5,2,m);
    wavefunctionc = @(t) spline(x,wavefunction1(n,:).*wavefunction1(n,:),t);
    wavefunction1(n,:) = wavefunction1(n,:)/sqrt(integral(wavefunctionc,x(1),x(length(x))));
    plot(x,wavefunction1(n,:));
    xlabel('x');
    ytext = strcat('$\psi_',num2str(n),'$');
    ylabel(ytext, 'Interpreter','Latex');
    titletext = strcat('the eigenstate of Infinite Square Well with sawtooth potential at level ',num2str(n));
    title(titletext);
end

%--------------------------------------------------------------------------
% 3. Infinite Square Well with triangle potential
%--------------------------------------------------------------------------
% doing the same method with the Infinite square well case but with adding
% the triangle potential to the k2 function handle

fprintf('3. Infinite Square Well Case with triangle potential hbar^2/m \n');
psi = zeros(1,n+1);
psi(1) = 0;
psi(2) = 1;
E = -1.;

% adding triangle potential to the k2 function handle and the rest of
% function handle.

k2 = @(x) 2.*(E-trianglepotential(1.5,x));
psi = numerov(psi(1),psi(2),x,k2,h);
psi = normalize(psi);
psiend = 0.0;
dE = 0.01;
tol = 1e-4;
wavefunction1 = zeros(5,length(psi));
Earray = zeros(1:10);
for n = 1:5
    while 1
        E = E + dE;
        k2 = @(x) 2.*(E-trianglepotential(1.5,x));
        psi1 = numerov(psi(1),psi(2),x,k2,h);
        psi1 = normalize(psi1);
        if psi(length(psi)).*psi1(length(psi1)) < 0
%             E = E - dE;
%             dE = dE/2.;
            % start from here
            E1 = E - dE;
            E2 = E;
             E1iter = E1;
                E2iter = E2;
            while abs(psi1(length(psi1))-psiend) > tol

                k2a = @(x) 2.*(E1iter - trianglepotential(1.5,x));
                k2b = @(x) 2.*(E2iter - trianglepotential(1.5,x));
                psia = numerov(0,1,x,k2a,h);
                psia = normalize(psia);
                psib = numerov(0,1,x,k2b,h);
                psib = normalize(psib);
                E3iter = E2iter - psib(length(psib)).*(E2iter - E1iter)/(psib(length(psib))-psia(length(psia)));
                k2c = @(x) 2.*(E3iter- trianglepotential(1.5,x));
                psi1 = numerov(0,1,x,k2c,h);
                psi1 = normalize(psi1);
                E2iter = E3iter;
                E1iter = E2iter;
            end
            psiresult = psi1;
            fprintf('The energy level at n = %d is %f hbar^2/m \n',n,E3iter);

            Earray(n) = E3iter;
            E = E2;
            % ends here
        end
        k2 = @(x) 2.*(E - trianglepotential(1.5,x));
        psi = numerov(0.0,1,x,k2,h);
        psi = normalize(psi);

        if abs(psi1(length(psi1))-psiend) < tol
            break;
        end

    end

    wavefunction1(n,:) = psiresult;
    dE = 0.01;
end

figure;
subplot(5,2,1:2:10);
hold on;
plot(x,trianglepotential(1.5,x),'LineWidth',4);
title('First 5 Energy Levels and the Infinite Square Well with triangle potential');
xlabel('x');
ylabel('Energy (hbar^2/m)');
for m=1:5
    plot(x,Earray(m).*ones(1,length(x)),'LineWidth',4);
end
legend('V_0','E_1','E_2','E_3','E_4','E_5');
hold off;
for n = 1:5
    m = 2.*n;
    subplot(5,2,m);
    wavefunctionc = @(t) spline(x,wavefunction1(n,:).*wavefunction1(n,:),t);
    wavefunction1(n,:) = wavefunction1(n,:)/sqrt(integral(wavefunctionc,x(1),x(length(x))));
    plot(x,wavefunction1(n,:));
    xlabel('x');
    ytext = strcat('$\psi_',num2str(n),'$');
    ylabel(ytext, 'Interpreter','Latex');
    titletext = strcat('the eigenstate of Infinite Square Well with triangle potential at level ',num2str(n));
    title(titletext);
end

