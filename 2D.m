%%%%%%%%%%%%%%%%%%%%%% Infinite potential well for 2D case %%%%%%%%%%%%%%%%%%%%
clc
close all
clear all

%%Since we are considering a particle ina square box by a potential
%%%energy V(x,y) ,V(x,y)=0 inside the box
%Inside the box entire energy is kinetic energy

%                   E = Ex+Ey = (px^2)/2m + (py^2)/2m
%
%Since the energy is a simple sum, solutions of the schrodinger equation
% for 2D case are siple products of solutions of the 1D case
% therefore wave functions have the form,
%                   psi(x,y)=A*sin(x(2m(Ex)/h^2)^0.5)*sin(y(2m(Ey)/h^2)^0.5) + 

%Using boundary conditions at psi(L,y)=psi(x,L)=0

L=10;   %width of the potential well
N= 1000; %Number of discrete points
x = linspace(0,L,N);
y = linspace(0,L,N);

n = 2;  %Number of waveforms
m = 2;  
%%%%%%%%%%%%%%%%Ploting the set of solution wave functions%%%%%%%%%%%%%%%%%
for i=1:n
    for j=1:m
        syms x y
        psi = (2/L)*sin((i*pi.*(x))/L).*sin((j*pi.*(y))/L);
        figure()
        title('\psi','fontSize',14);
        ezsurf(psi,[0,L,0,L])
        grid on
    end
end
%%%%%%%%%%%%%%%%Probability Density functions for each solution%%%%%%%%%%%
for i=1:n
    for j=1:m
        syms x y
        Pr = abs(((2/L)*sin((i*pi.*(x))/L).*sin((j*pi.*(y))/L))^2);
        figure()
        title('Pr(x,y)','fontSize',14);
        ezsurf(Pr,[0,L,0,L])
        grid on
    end
end