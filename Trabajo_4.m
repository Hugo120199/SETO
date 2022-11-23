%% Trabajo 4
    % A partir de la pala y el report de NREL 5MW
clc
clear all
close all


%% Datos   
    % Condiciones de la pala
Psi = pi/2;
Phi = pi/2;
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
    [~,~,~,~,~,cl] = ff('cl');
    [~,~,~,~,~,cm] = ff('cm');
    for j = 1:length(cd)
        if abs(alpha2(1,i)) <= abs(alpha(j,1))*1.05 & abs(alpha2(1,i)) >= abs(alpha(j,1))*0.95 
            cd_x(i) = cd(j,1);
            break
        else
            cd_x(i) = 0;
        end
    end
    for j = 1:length(cl)
        if abs(alpha2(1,i)) <= abs(alpha(j,1))*1.05 & abs(alpha2(1,i)) >= abs(alpha(j,1))*0.95 
            cl_x(i) = cl(j,1);
            break
        else
            cl_x(i) = 0;
        end
    end
    for j = 1:length(cm)
        if abs(alpha2(1,i)) <= abs(alpha(j,1))*1.05 & abs(alpha2(1,i)) >= abs(alpha(j,1))*0.95 
            cm_x(i) = cd(j,1);
            break
        else
            cm_x(i) = 0;
        end
    end

end

%% Cálculo Ca(r) y Fb,zh1(r)
    % Para pala sin rotación (Omega = 0 rad/s) la expresión de Ca(r) es: 
    % Ca(r) rho*U*c(r)*cd(r)

for i = 1:length(r)
    Ca(i) = rho*U1*c(i)*cd_x(i);
    fb_zh1(i) = 0.5*U1*Ca(i);
end
figure(1)
subplot(1,2,1)
plot(r,Ca,'b-')
grid on;
xlabel('$r$ [m]','interpreter','latex','fontsize',12)
ylabel('$C_a (r)$ [Ns/m^{2}]','interpreter','latex','fontsize',12)

subplot(1,2,2)
plot(r,fb_zh1,'r-')
grid on;
xlabel('$r$ [m]','interpreter','latex','fontsize',12)
ylabel('$F_{b,zH1} (r)$ [N/m]','interpreter','latex','fontsize',12)

%% Análisis modal
    % Características estructurales de la pala
pala = beamNREL_5MW49tp(@bladeNREL_5MW49tp,thetaC);
m = pala.mass;

    % Solución modal
[omega1,w,dw] = solveEigen(pala,1,0);
psi1 = w;
    % Vector deformación del primer modo
figure(2)
plot(r,w)
grid on;
xlabel('$r$ [m]','interpreter','latex','fontsize',12)
ylabel('$w (r)$ [m]','interpreter','latex','fontsize',12)

    % Cálculo de las características modales
m1 = trapz(r,psi1.^2.*m);
ca1 = trapz(r,psi1.^2.*Ca');
chi_a1 = ca1/(2*m1*omega1/(2*pi));


%% 
function [Ay,Az,Ayz] = pAxes2aAxes(A1,A2,alpha)
    % According to Bauchau page 234 equations (6.34a,b,c) where for the
    % definition of principal axes we have set A12 = 0
    
    Ay       = 0.5*(A1 + A2) + 0.5*(A1 - A2).*cos(2*alpha);
    Az       = 0.5*(A1 + A2) - 0.5*(A1 - A2).*cos(2*alpha);
    Ayz      = 0.5*(A1 - A2).*sin(2.*alpha);
end