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
xlabel('$r [m]$','interpreter','latex','fontsize',12)
ylabel('$C_a (r)$ [Ns/m^{2}]$','interpreter','latex','fontsize',12)

subplot(1,2,2)
plot(r,fb_zh1,'r-')
grid on;
xlabel('$r [m]$','interpreter','latex','fontsize',12)
ylabel('$F_{b,zH1} (r) [N/m]$','interpreter','latex','fontsize',12)

%% Análisis modal
    % Características estructurales de la pala
pala = beamNREL_5MW49tp(@bladeNREL_5MW49tp,thetaC);
m = pala.mass;

    % Solución modal
[omega1,w,dw] = solveEigen(pala,1,0);
psi1 = w;
f1 = omega1/(2*pi);
    % Vector deformación del primer modo
figure(2)
plot(r,w)
grid on;
xlabel('$r$ [m]','interpreter','latex','fontsize',12)
ylabel('$w (r)$ [m]','interpreter','latex','fontsize',12)

    % Cálculo de las características modales
m1 = trapz(r,psi1.^2.*m);
ca1 = trapz(r,psi1.^2.*Ca');
chi_a1 = ca1/(2*m1*(f1)); %esto XD no da muy bien pero el resto cuadra
chi_s1 =  0.004775; % Considerar el amortiguamiento estructural proporcionado en la referencia Jonkman et al. (2009).
cs1 = 2*m1*chi_s1*omega1;
c1 = ca1 + cs1;
chi_1 = chi_a1 + chi_s1;

    % Imprimimos los resultados por pantalla:
disp("----------------------------------------------------------------------");
disp("f_1 [Hz] | m_1 [kg] | c_1 [kg/s] | chi_1a [%] | chi_1s [%] | chi_1 [%]");
disp(string(omega1/(2*pi)) + " | " + string(m1) + " | " + string(c1) + " | " + string(chi_a1*100) + " | " + string(chi_s1*100) + " | " + string(chi_1*100));

    % Función de transferencia del libro que nos piden el bode:

H1 = tf(1,m1*[1, 2*chi_1*omega1, omega1^2]);
[mag,phase,wout] = bode(H1);

figure()
subplot(2,1,1)
semilogx(wout/(2*pi),20*log10(squeeze(mag(1,1,:))),"-","Color",'black');
ylabel("$Magnitud [dB]$",'Interpreter','latex');
xlim([wout(1)/(2*pi) wout(end)/(2*pi)]);
grid minor;

subplot(2,1,2)
semilogx(wout/(2*pi),squeeze(phase(1,1,:)),"-","Color",'black');
ylabel("$Fase [^{\circ}]$",'Interpreter','latex');
xlabel("$f [Hz]$",'Interpreter','latex');
xlim([wout(1)/(2*pi) wout(end)/(2*pi)]);
grid minor;

    % Comparación de f1 con los distintos modelos (esta en el libro esto también)

 f1_fast = 1.0793; error_fast = abs(f1_fast - f1)/f1;
 f1_adams = 1.074; error_adams = abs(f1_adams - f1)/f1;

 disp("---------------------------------------");
 disp("           |   f1   | f1_fast | f1_adams ");
 disp("Valor [Hz] | " + string(f1) + " | " + string(f1_fast) + " | " + string(f1_adams));
 disp("Error [%]  | " + " ----- | " + string(error_fast*100) + " | " + string(error_adams*100));

 %% Distribución promedio de desplazamiento y momento flector de batimiento

%% 
function [Ay,Az,Ayz] = pAxes2aAxes(A1,A2,alpha)
    % According to Bauchau page 234 equations (6.34a,b,c) where for the
    % definition of principal axes we have set A12 = 0
    
    Ay       = 0.5*(A1 + A2) + 0.5*(A1 - A2).*cos(2*alpha);
    Az       = 0.5*(A1 + A2) - 0.5*(A1 - A2).*cos(2*alpha);
    Ayz      = 0.5*(A1 - A2).*sin(2.*alpha);
end