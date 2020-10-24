%This is a program coding to solve Triangular fuzzy problem
%Depending upon the problem the matrices MatA and U0 will vary
clear all
clc
format long
m=1000;
t=1.0;
alp=[0:0.1:1];
bet=[0:0.1:1];
ShowU=zeros(8,length(alp));
for k=1:length(alp)
X0=[26 28 29 27 28 30]';
Y0=[15 17 20 14 17 19]';
A=3;
B=4;
I=eye(8);
MatA=[0 0 0 0 A 0 0 0;
      0 0 0 0 0 A 0 0;
      0 0 0 0 0 0 A 0;
      0 0 0 0 0 0 0 A;
      B 0 0 0 0 0 0 0;
      0 B 0 0 0 0 0 0;
      0 0 B 0 0 0 0 0;
      0 0 0 B 0 0 0 0];
S=MatA;
%Inti=[0 0 0 0 0 0 0 0]'
%initial conditions definitions
U0=[X0(1)+alp(k)*(X0(2)-X0(1)); 
      X0(3)-alp(k)*(X0(3)-X0(2));
      X0(2)-bet(k)*(X0(2)-X0(4));
      X0(2)+bet(k)*(X0(6)-X0(2));
      Y0(1)+alp(k)*(Y0(2)-Y0(1));
      Y0(3)-alp(k)*(Y0(3)-Y0(2));
      Y0(2)-bet(k)*(Y0(2)-Y0(4));
      Y0(2)+bet(k)*(Y0(6)-Y0(2))];
U=zeros(8,t*m);
G=(I-1/(2*m).*S)\(1/m.*MatA*U0);
U(:,1)=G+U0;
%%%%%%%%%%%%%
for i=2:t*m
G=(I-1/(2*m).*S)\(1/m.*MatA*U(:,i-1));
U(:,i)=G+U(:,i-1);
end
U(:,t*m);
ShowU(:,k)=U(:,t*m);
end
ShowX=ShowU(1:4,:);
ShowY=ShowU(5:8,:);
FinalX=[alp' ShowX']
FinalY=[alp' ShowY']
%plot(ShowX,alp')
%plot(ShowY,alp')
figure;plot(ShowX,alp)
figure;plot(ShowY,alp)