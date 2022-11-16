%Buscar NREL-5MW --> Amortiguamiento estructural

Blade_geo = bladeNREL_5MW49tp.m; %Características geométricas
Blade_st = beamNREL_5MW49tp.m; %Características estructurales

%%Datos
Psi = pi/2;
U1 = 12;
rho = 1.225;
Iu = 0.15;
Lu = 350;
thetaC = pi/2;

function [Ay,Az,Ayz] = pAxes2aAxes(A1,A2,alpha)
% According to Bauchau page 234 equations (6.34a,b,c) where for the
% definition of principal axes we have set A12 = 0

Ay       = 0.5*(A1 + A2) + 0.5*(A1 - A2).*cos(2*alpha);
Az       = 0.5*(A1 + A2) - 0.5*(A1 - A2).*cos(2*alpha);
Ayz      = 0.5*(A1 - A2).*sin(2.*alpha);