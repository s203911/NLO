clear all;
close all;
clc;

% defining basi constats
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
a1=sqrt(P1/(hbar*omega1*pi*W0^2));
a2=sqrt(P2/(hbar*omega2*pi*W0^2));
a3=sqrt(P3/(hbar*omega3*pi*W0^2));

astart=[a1, a2, a3];
[Z,A] = ode45(@(z,a) secondorder(z,a,K,g,dk),[0 1],astart);

Pgen1=abs(A(:,1)).^2*hbar*omega1*pi*W0^2;
Pgen2=abs(A(:,2)).^2*hbar*omega2*pi*W0^2;
Pgen3=abs(A(:,3)).^2*hbar*omega3*pi*W0^2;

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



