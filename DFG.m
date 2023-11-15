close all
clear all
clc;
% definition of units
global m nm kg sec A K mol cd J eV mass_e h_bar epsilon0 c pm V
m = 1;
nm = 1E-9 * m;
pm = 1E-12 * m;
kg = 1;
sec = 1;
A = 1;
K = 1;
V = 1; 
mol = 1;
cd = 1; % candela
J = 1;
eV = 1.60218E-19 * J;
mass_e = 9.10938E-31 * kg;
h_bar = (6.62607E-34) / (2 * pi) * (J * sec);
epsilon0 = 8.8541878 * A^2 * sec^4 / (kg * m^3);
c = 299792458 * m / sec;


%% Constants from setup
% the coupling w's and lambda's
lambda3 = 1064 * nm; % pump 
lambda1 = 1550 * nm; % signal 
lambda2 = 1/( 1/lambda3 - 1/lambda1 ); % idler, generated
omega1=(c/lambda1)*2*pi;
omega2=(c/lambda2)*2*pi;
omega3=(c/lambda3)*2*pi;



% The refractive index based on the sellmeier equations
T = 297.66 * K; % room T at 24.5 C
lambda = 1.55 * 10^(-6) * m; % no??
n = neo(lambda,T); % refractive index for extraordinary refractive index

% from this, we can compute eta as
eta = 377/n;

chi_eff= 17.2 * pm / V; % effective chi^(2) for DFG for LnB

K=1;
dk=0;

% computing g_dfg^2=phi_w3(0)*g^2
g = epsilon0 * chi_eff * sqrt((1/2) * eta^3 *h_bar * omega1 * omega2 * omega3)

W0=100e-6; % power 





%%


P1=100;
P2=1;
P3=1;
% CHANGE NAME OF a1 a2 a3
a1=sqrt(P1/(hbar*omega1*pi*W0^2));
a2=sqrt(P2/(hbar*omega2*pi*W0^2));
a3=sqrt(P3/(hbar*omega3*pi*W0^2));

astart=[a1, a2, a3];
[Z,A] = ode45(@(z,a) secondorder(z,a,K,g,dk),[0 1],astart);

Pgen1=abs(A(:,1)).^2*h_bar*omega1*pi*W0^2;
Pgen2=abs(A(:,2)).^2*h_bar*omega2*pi*W0^2;
Pgen3=abs(A(:,3)).^2*h_bar*omega3*pi*W0^2;

plot(Z,Pgen1-100)
hold on
plot(Z,Pgen2)
hold on
plot(Z,Pgen3)
hold on
plot(Z,Pgen1+Pgen2+Pgen3-100)
legend('P1-100','P2','P3','Sum-100')

%% Opgave b
eps0=8.85e-12;
hbar=1.054e-34;
xeff=5e-12;
c=3e8;
K=1;
dk=0;
n=2.2;
eta=377/n;
lambda1=1342e-9;
lambda2=1064e-9;
lambda3inv=(1/lambda1+1/lambda2);
lambda3=1/lambda3inv;
omega1=(c/lambda1)*2*pi;
omega2=(c/lambda2)*2*pi;
omega3=(c/lambda3)*2*pi;
g=eps0*K*xeff*sqrt(1/2*eta^3*hbar*omega1*omega2*omega3);
W0=100e-6;


P1=100;
P2=1;
P3=1;

for q=[1 2 3 4]
    a1=sqrt(P1/(hbar*omega1*pi*W0^2));
    a2=sqrt(P2/(hbar*omega2*pi*W0^2));
    a3=sqrt(P3/(hbar*omega3*pi*W0^2))*exp(1i*q*pi/2);
    
    astart=[a1, a2, a3];
    [Z,A] = ode45(@(z,a) secondorder(z,a,K,g,dk),[0 1],astart);

    Pgen1=abs(A(:,1)).^2*hbar*omega1*pi*W0^2;
    Pgen2=abs(A(:,2)).^2*hbar*omega2*pi*W0^2;
    Pgen3=abs(A(:,3)).^2*hbar*omega3*pi*W0^2;

    subplot(2,2,q)
    plot(Z,Pgen1-100)
    hold on
    plot(Z,Pgen2)
    hold on
    plot(Z,Pgen3)
    hold on
    plot(Z,Pgen1+Pgen2+Pgen3-100)
    legend('P1-100','P2','P3','Sum-100')
end

%% Opgave c
lambda1um = 1.342;
lambda2um = 1.064;
lambda3invum=(1/lambda1um+1/lambda2um);
lambda3um=1/lambda3invum;

ne1 = 4.5820+(0.099169/((lambda1um)^2-0.04443))-0.02195*(lambda1um)^2;
ne2 = 4.5820+(0.099169/((lambda2um)^2-0.04443))-0.02195*(lambda2um)^2;
ne3 = 4.5820+(0.099169/((lambda3um)^2-0.04443))-0.02195*(lambda3um)^2;

phasemismatch = ne1/lambda1um+ne2/lambda2um-ne3/lambda3um

QPMlength = 2*pi/(2*phasemismatch)

%% Opgave d

xeff = 5e-12;
g=eps0*K*xeff*sqrt(1/2*eta^3*hbar*omega1*omega2*omega3);

P1=100;
P2=1;
P3=0;
a1=sqrt(P1/(hbar*omega1*pi*W0^2));
a2=sqrt(P2/(hbar*omega2*pi*W0^2));
a3=sqrt(P3/(hbar*omega3*pi*W0^2));

astart=[a1, a2, a3];
[Z,A] = ode45(@(z,a) secondorder(z,a,K,g,dk),[0 1],astart);

Pgen3=abs(A(:,3)).^2*hbar*omega3*pi*W0^2;

plot(Z,Pgen3)
hold on

xeff = 25e-12;
g=eps0*K*xeff*sqrt(1/2*eta^3*hbar*omega1*omega2*omega3);

astart=[a1, a2, a3];
[Z,A] = ode45(@(z,a) secondorder(z,a,K,g,dk),[0 1],astart);

Pgen3=abs(A(:,3)).^2*hbar*omega3*pi*W0^2;

plot(Z,Pgen3)



