%% Trabajo 4
    % A partir de la pala y el report de NREL 5MW
clc
clear all
close all


%% Datos   
    % Condiciones de la pala
Psi = pi/2;
Phi = 0;
U1 = 12;
rho = 1.225;
Iu = 0.15;
Lu = 350;
thetaC = pi/2;      % Pala en bandera
R = 63;

    % Pala (características) 
[nn, ndc] = bladeNREL_5MW49tp('cadvec');
[nn, x] = bladeNREL_5MW49tp('xvec');
c = ndc.*R;
r = x.*R;
[nn, c_aero] = bladeNREL_5MW49tp('airfoilvec');
[nn, thetaG] = bladeNREL_5MW49tp('twistvec');


for i = 1:length(c_aero)
    alpha2(i) = Phi-thetaC-thetaG(i);
    ff = c_aero{i}; 
    [~,~,~,~,alpha,cd] = ff('cd');
    for j = 1:length(cd)
        if alpha2(i) <= alpha(j)*1.05 && alpha2(i) >= alpha(j)*0.95 
            cd_x(i) = cd(j);
        else
            cd_x(i) = 0;
        end
    end


end

%     [alpha, cd_x] = ff('cd');
%     [alpha, cm_x] = ff('cm');
%     cl(i) = cl_x;
%     cd(i) = cd_x;
%     cm(i) = cm_x;
% end

%% Cálculo Ca(r) y Fb,zh1(r)
    % Para pala sin rotación (Omega = 0 rad/s) la expresión de Ca(r) es: 
    % Ca(r) rho*U*c(r)*cd(r)

for i = 1:length(r)
    Ca(i) = rho*U1*c(i)*cd(i);
    fb_zh1(i) = 0.5*U1*Ca(i);
end
figure(1)
plot(r,Ca,'b-',r,fb_zh1,'r-')
xlabel('$r$ [m]','interpreter','latex','fontsize',12)
%ylabel('$C_a(r)$ $[Ns/m^2]$','interpreter','latex','fontsize',12)
legend('Ca(r) [Ns/m^{2}]','$F_{b,zH1}$ [N/m]','interpreter','latex')

%% Análisis modal



%% 
function [Ay,Az,Ayz] = pAxes2aAxes(A1,A2,alpha)
    % According to Bauchau page 234 equations (6.34a,b,c) where for the
    % definition of principal axes we have set A12 = 0
    
    Ay       = 0.5*(A1 + A2) + 0.5*(A1 - A2).*cos(2*alpha);
    Az       = 0.5*(A1 + A2) - 0.5*(A1 - A2).*cos(2*alpha);
    Ayz      = 0.5*(A1 - A2).*sin(2.*alpha);
end