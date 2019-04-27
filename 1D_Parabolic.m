%Chandula nethmal
%jan2019

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%% Defining parameters and constants %%%%% %%%%%%%%%%%%%
E1=[];
h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);                  
e=1.602176487E-19;              %% electron charge [C]
me=9.10938188E-31;              %% electron mass [kg]

dz=2E-10;               % resolution of the grid [m]
n=4;                   % number of solution asked 
Mass = 0.067;           % effective mass, constant over all the structure...
C=0.1;                % scaling factor to plot the wave function [Without Dimension]
a=1e15;
L=50e-9;
z=0:dz:L;
%%%%%%%%%%%%%%%%%%%%%%%%% Potential definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vb=0.3;                 % potential barrier [eV]
Nz=length(z);
V0=zeros(1,Nz);
V0= a*(z-L/2).^2;
V0(V0>Vb)=Vb;
V0=V0-min(V0);

figure(1)
plot(z*1e9,V0, 'b-','linewidth',2)
xlabel('z (nm)');
ylabel('V(z)');
ylim([min(V0)-0.05 max(V0)+0.1]);
title('Parabolic potential');

%%%%%%%%%%%%%%%%%%%%% Building the operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DZ2 =(-2)*diag(ones(1,Nz),0) + (1)*diag(ones(1,Nz-1),-1) + (1)*diag(ones(1,Nz-1),1);
DZ2=DZ2/dz^2;

%%%%%%%%%%%%%%%%%%%% Building of the Hamiltonian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = -hbar^2/(2*me*Mass)*DZ2  +  diag(V0*e);
H = sparse(H);   % compact the H ingnoring non zeros
[psi,Energy] = eigs(H,n,'SM');  % Directly finding eigen values of the hamiltonian
E = diag(Energy)/e ;
E=real(E);
E=E(end:-1:1);
E1=E;
E=zeros(n,4);
E(1:length(E1),1)=E1;

%%%%%%%%%%%%%%%%%% Normalization of the Wavefunction %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
for i=1:n
    psi(:,i)=psi(:,i)/sqrt(trapz(z',abs(psi(:,i)).^2));  % normalisation at 1
end

psi=psi(:,end:-1:1);
psi1=psi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(strcat('E(eV)='))
display(strcat(num2str(E)))

figure(2)
xlabel('z (nm)');
ylabel('\psi ');
hold on


for i=1:length(E1)
  subplot(4,1,i);
  plot(z*1e9,psi1(:,i).*1e-4,'r-','linewidth',2)
  xlabel('z (nm)','fontSize',10);
  ylabel('\psi','fontSize',14);
end

for i=1:length(E1)
  psi1(:,i)=abs(psi1(:,i)).^2/max(abs(psi1(:,i)).^2)*C + E1(i) % normalisation for the plotting
end

figure(3)
subplot(1,1,1,'fontsize',15);
hold on;
%grid on;
plot(z*1e9,V0, 'b-','linewidth',2)

xlabel('z (nm)');
ylabel('Energy (eV)');
ylim([min(V0)-0.05 max(V0)+0.1]);
title('Probability Functions');
for i=1:length(E1)
  plot(z*1e9,psi1(:,i),'r-','linewidth',2)
end
