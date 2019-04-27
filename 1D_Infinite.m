%%%%%%%%%%%%%%%%%%%%%% Infinite potential well for 1D case %%%%%%%%%%%%%%%%%%%%
clc
close all
clear all

%%Since we are considering a particle in a infinite potential well
%%% V(z)=0 inside the well
%Inside the well entire energy is kinetic energy

% E = (pz^2)/2m
% therefore wave functions have the form
%psi(z)=A*sin(x(2m(Ex)/h^2)^0.5) 

%Using boundary conditions at psi(L,y)=psi(x,L)=0

L=10;   %width of the potential well in nm
N=1000;
z = linspace(0,L,N);
%%%%%%%%%%%%%%%%%%%%Defining the solution%%%%%%%%%%%%%%%%%%
n = 4;  %Number of waveforms
figure
for i=1:n
    
    psi = (2/L)*sin((i*pi.*(z))/L);
    subplot(4,1,i) 
    plot(z,psi,'b','linewidth',2);
    xlabel('z (nm)','fontSize',14);
    ylabel('\psi','fontSize',14);
end

figure
%%%%%%%%%%%%%%%Plotting the probability Distribution%%%%%%%%
for i=1:n
    
    P = abs(((2/L)*sin((i*pi.*(z))/L)).^2);
    subplot(4,1,i) 
    plot(z,P,'b','linewidth',2);
    xlabel('z (nm)','fontSize',14);
    ylabel('P(z)','fontSize',14);
end
