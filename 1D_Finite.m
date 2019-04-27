%%%%%%%%%%%%%% Schrodinger Eqaution for 1D finite potential well %%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
%Chandula nethmal
%jan2019
clc
close all
clear all

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
q=1.602176487e-19;              %% electron charge [C]
m=9.10938188E-31;

Lz=1e-9; %9nm
N=200;  %200 Samples
Vo=1;   %1eV
V=zeros(1,N+1);
z=linspace(-2*Lz,2*Lz,N+1);
%%%%%%%%%%%%%%%%Definig the potential function%%%%%%%%%%%%%%
V(1,1:N/4+1)=Vo;
V(1,3*(N/4):N+1) = Vo;
figure(1)
%subplot(2,1,l)
plot(z,V,'linewidth',2)
axis([-2*Lz 2*Lz 0 Vo+0.2])
xlabel('z (m)','fontSize',14);
ylabel('V (eV)','fontSize',14);
%-(h^2/2m)d2Psi/dz2 + V(z).Psi=E.Psi
%Applying Schrodinger's equation to diffrent regions of the well

%Psi =G*exp(K*z)                 , -100<z<-50
%Psi =A*sin(k*z) + B*cos(k*z	 , -50 <z<50
%Psi =F*exp(K*z)                 ,  50 <z<100    

%Using Boundary conditions for the above three equations we can simplify the system
%According to the nature of A,B,G and F coefficients, we will have symmetric and antisymmetric solutions 
%tan(k*2Lz/2)=K/k if F=G
%-cot(k*2Lz/2)=K/k if F=-G
Vo=Vo*q;
Einf=((hbar^2)*pi*pi)/(2*m*(2*Lz)^2) % First Energy level of the infinite well
vo = Vo/Einf
Epsilon =linspace(0,20,200);%0:0.01:15;

%For the Symmetric solutions we ahve
%(Epsilon^0.5)*tan((pi/2)*(Epsilon^0.5)) = (vo -Epsilon)^0.5 

%For the Antiymmetric solutions we ahve
%-(Epsilon^0.5)*cot((pi/2)*(Epsilon^0.5)) = (vo -Epsilon)^0.5

% For each eqation We need to find the intersections left side curves with the right side curve
E1=(vo -Epsilon).^0.5;
E2=(Epsilon.^0.5).*tan((pi/2)*(Epsilon.^0.5));
E3=-(Epsilon.^0.5).*cot((pi/2)*(Epsilon.^0.5));
E4=Epsilon;

figure(2)
plot(Epsilon,E1,'r','linewidth',2);
hold
plot(Epsilon,E2,'b','linewidth',2);
plot(Epsilon,E3,'g','linewidth',2);
axis([0 20 0 10]);
xlabel('Epsilon','fontSize',14);
ylabel('E1-R,E2-G,E3-B','fontSize',14);
grid on
%Here we find that there are four points in the plot are intersecting
%Two of them are symmetric and other two are antisymmetric

even_index_intersection= find(abs(E1-E2)<=0.1);
odd_index_intersection= find(abs(E1-E3) <=0.18);

even_x_val=Epsilon(even_index_intersection)
odd_x_val=Epsilon(odd_index_intersection)
xval=sort([even_x_val,odd_x_val]);
a=size(xval,2);
xval=xval(1:a-1)
%define coefficients for wave functions
A=1;
psi=zeros(1,201);
figure(3)
z=z/(2*Lz);
for i=1:a-1;
    e0=xval(i);
    Cl= cos(pi*(e0^0.5)*0.5)/exp(-pi*((vo-e0)^0.5)*0.5);
    Sl= sin(pi*(e0^0.5)*0.5)/exp(-pi*((vo-e0)^0.5)*0.5);
    if mod(i+1,2)==0
        psi(1:N/4+1)=A*Cl*exp(pi*((vo-e0)^0.5).*z(1:(N/4)+1));
        psi(N/4+1:3*(N/4)+1)=A*cos(pi*((e0)^0.5).*z(N/4+1:3*(N/4)+1));
        psi(3*(N/4)+1:N+1)=A*Cl*exp(-pi*((vo-e0)^0.5).*z(3*(N/4)+1:N+1));
    else
        psi(1:N/4+1)=-A*Sl*exp(pi*((vo-e0)^0.5).*z(1:(N/4)+1));
        psi(N/4+1:3*(N/4)+1)=A*sin(pi*((e0)^0.5).*z(N/4+1:3*(N/4)+1));
        psi(3*(N/4)+1:N+1)=A*Sl*exp(-pi*((vo-e0)^0.5).*z(3*(N/4)+1:N+1));
    end
    subplot(4,1,i)
    %plot(z,V,'linewidth',2)
    %hold on
    plot(z,psi,'linewidth',2)
    xlabel('z (nm)','fontSize',14);
    ylabel('\psi','fontSize',14);
end
