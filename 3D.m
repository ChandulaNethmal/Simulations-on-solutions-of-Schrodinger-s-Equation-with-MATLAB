%%%%%%%%%%%%%%%%%%%%%% particle ina infintely deep three dimentional quantum well %%%%%%%%%%%%%%%%%%%%
clc
close all
clear all

%%Since we are considering a particle in a 3D well by a potential
%%%V(x,y,z) ,
% V(x,y,z)=0;           inside the box,
% V(x,y,z)= inifinity;  otherwise

%This function can be written as a sum of three 1D functions

%                       V(r) = V(x) + V(y) + V(z)

% Since the potential is a simple sum, solutions of the schrodinger equation
% for 3D case are simple products of solutions of the 1D case
% therefore wave functions have the form
%                       psi(r) = psi(x)*psi(y)*psi(z)

%  Energies also will be added up together to find Energy for 3D case
%                       psi(r)=psi(x,y,z)=A*sin(x(2m(Ex)/h^2)^0.5)*sin(y(2m(Ey)/h^2)^0.5)*sin(z(2m(Ez)/h^2)^0.5) 
%            
%Using boundary conditions at psi(L,y)=psi(x,L)=0

L=10;   %Width of the potential well
dx=0.5; dy=0.5; dz=0.5;
% define an array of values for x y and z
x = 0:dx:L; 
y = 0:dy:L; 
z = 0:dz:L; 
% creating a meshgrid
[Xmesh,Ymesh,Zmesh] = meshgrid(x,y,z);
% compute the function value at the meshgrid points
m=2; n=2;p=2;
for i=1:m
    for j=1:n
        for k=1:p
            f = ((8/(L^3))^0.5)*sin(i*pi*Xmesh/L).*sin(j*pi*Ymesh/L).*sin(k*pi*Zmesh/L);
            figure
            xslice=[2,5,8]; yslice=[3,6]; zslice = [5];     %Defining slicing planes for 3 axis
            slice(Xmesh,Ymesh,Zmesh,f,xslice,yslice,zslice)
            title(sprintf('psi(x,y,z), when nx,ny,nz= %d, %d, %d',i,j,k));
            xlabel('x-axis');
            ylabel('y-axis');
            zlabel('z-axis');
            colormap jet
            colorbar
        end
    end
end