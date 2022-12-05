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
    cd_x(i) = interp1(alpha,cd,alpha2(i));
    cl_x(i) = interp1(alpha,cd,alpha2(i));
    cm_x(i) = interp1(alpha,cd,alpha2(i));

    % Gran hazaña, pero hemos visto que había soluciones mejores
%     for j = 1:length(cd)
%         if abs(alpha2(1,i)) <= abs(alpha(j,1))*1.05 & abs(alpha2(1,i)) >= abs(alpha(j,1))*0.95 
%             cd_x(i) = cd(j,1);
%             break
%         else
%             cd_x(i) = 0;
%         end
%     end
%     for j = 1:length(cl)
%         if abs(alpha2(1,i)) <= abs(alpha(j,1))*1.05 & abs(alpha2(1,i)) >= abs(alpha(j,1))*0.95 
%             cl_x(i) = cl(j,1);
%             break
%         else
%             cl_x(i) = 0;
%         end
%     end
%     for j = 1:length(cm)
%         if abs(alpha2(1,i)) <= abs(alpha(j,1))*1.05 & abs(alpha2(1,i)) >= abs(alpha(j,1))*0.95 
%             cm_x(i) = cd(j,1);
%             break
%         else
%             cm_x(i) = 0;
%         end
%     end

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
xlabel('$r$ [m]','interpreter','latex','fontsize',14)
ylabel('$C_a$ $(r)$ [Ns/m$^{2}$]','interpreter','latex','fontsize',14)

subplot(1,2,2)
plot(r,fb_zh1,'r-')
grid on;
xlabel('$r$ [m]','interpreter','latex','fontsize',14)
ylabel('$F_{b,zH1}$ $(r)$ [N/m]','interpreter','latex','fontsize',14)

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
xlabel('$r$ [m]','interpreter','latex','fontsize',14)
ylabel('$\psi_1$','interpreter','latex','fontsize',14)

    % Cálculo de las características modales
m1 = trapz(r,psi1.^2.*m);
ca1 = trapz(r,psi1.^2.*Ca');
chi_a1 = ca1/(2*m1*(f1));   % esto XD no da muy bien pero el resto cuadra
chi_s1 =  0.004775;     % Considerar el amortiguamiento estructural proporcionado en la referencia Jonkman et al. (2009).
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

figure(3)
subplot(2,1,1)
semilogx(wout/(2*pi),20*log10(squeeze(mag(1,1,:))),"-","Color",'black');
ylabel("$Magnitud$ [dB]",'Interpreter','latex','fontsize',14);
xlabel("$f$ [Hz]",'Interpreter','latex','fontsize',14);
xlim([wout(1)/(2*pi) wout(end)/(2*pi)]);
grid minor;

subplot(2,1,2)
semilogx(wout/(2*pi),squeeze(phase(1,1,:)),"-","Color",'black');
ylabel("$Fase$ $[^{\circ}]$",'Interpreter','latex','fontsize',14);
xlabel("$f$ [Hz]",'Interpreter','latex','fontsize',14);
xlim([wout(1)/(2*pi) wout(end)/(2*pi)]);
grid minor;

    % Comparación de f1 con los distintos modelos (esta en el libro esto también)
f1_fast = 1.0793;  %f1_fast = 0.6664; ??? % Esto es por estar en bandera
error_fast = abs(f1_fast - f1)/f1;
f1_adams = 1.074;  %f1_adams = 0.6292; ???
error_adams = abs(f1_adams - f1)/f1;

disp("---------------------------------------");
disp("           |   f1   | f1_fast | f1_adams ");
disp("Valor [Hz] | " + string(f1) + " | " + string(f1_fast) + " | " + string(f1_adams));
disp("Error [%]  | " + " ----- | " + string(error_fast*100) + " | " + string(error_adams*100));

%% Distribución promedio de desplazamiento y momento flector de batimiento
    % Distribuación promedio de W 
Fb_zh1_prom = 1./(m1*(omega1)^2).*trapz(r,fb_zh1'.*psi1);
W_prom = psi1.*Fb_zh1_prom;
    % Pequeños--> empotrada (pequeños desplazamientos pero mayores momentos)I 
    % y fuerzas asociadas a la posición de bandera (no gira y le entra
    % aire alineado a los perfiles---> resistencia pequeña)
figure(4)
plot(r,W_prom,'k-')
grid on;
xlabel('$r [m]$','interpreter','latex','fontsize',14)
ylabel('$W (r)$ [m]','interpreter','latex','fontsize',14)

    % Distribución promedio de momento flector (esto está con pinzas)
for i = 1:length(r)
    Mb_yh1_prom(i) = 1./(m1*(omega1)^2).*trapz(r,fb_zh1'.*psi1.*r(i));
end
figure(5)
plot(r,Mb_yh1_prom,'k-')
grid on;
xlabel('$r [m]$','interpreter','latex','fontsize',14)
ylabel('$Mb_{yH1} (r)$ [Nm]','interpreter','latex','fontsize',14)

%% Varianza de desplazamiento y de momento flector de batimiento

    % Para la condición de viento S0C1, representar la dependencia funcional de la varianza
% del desplazamiento en la punta, σ2 w,t, con la intensidad de turbulencia, Iu, y con la
% longitud de escala, Lu. Comentar y discutir dicha dependencia funcional.


Lu = linspace(50,500,50);
Iu = linspace(0.1,0.2,5);

for l=1:1:49
    v(l)=psi1(l)*Ca(l);
    v_2(l) = m(l)*psi1(l)*r(l);
end

for i = 1:length(Lu)
    for j = 1:length(Iu)
        
        k_co = trapz(r,v);
        k_mb1 = trapz(r,v_2);
        So = 4*(Iu(j)^2)*Lu(i)*U1*(2/3)^(5/3);
        sigma_u = (U1*Iu(j))^2;
        sigma_wt(i,j) = (k_co.^2)*((sigma_u/(16*pi*(m1.^2)*f1.^4)) + So/(64*pi*chi_1*(m1.^2)*f1.^3));
        sigma_mbr(i,j) = ((2*pi*f1)^2)*(k_mb1^2)*(sigma_wt(i,j));
    end
end

LW = 1.5; FS = 18;

figure()
plot(Lu,sigma_wt(:,1),'-','LineWidth',LW); hold on
plot(Lu,sigma_wt(:,2),'-','LineWidth',LW); hold on
plot(Lu,sigma_wt(:,3),'-','LineWidth',LW); hold on
plot(Lu,sigma_wt(:,4),'-','LineWidth',LW); hold on
plot(Lu,sigma_wt(:,5),'-','LineWidth',LW); hold on

xlabel("$Lu [m]$",'Interpreter','latex','FontSize',FS);
ylabel("$\sigma_{w,t}^2 [m^2]$",'Interpreter','latex','FontSize',FS);
legend("$Iu = " + string(Iu) + " $",'Interpreter','latex','Location','northwest','FontSize',FS);
grid on;

figure()
plot(Iu,sigma_wt(1,:),'-','LineWidth',LW); hold on
plot(Iu,sigma_wt(10,:),'-','LineWidth',LW); hold on
plot(Iu,sigma_wt(25,:),'-','LineWidth',LW); hold on
plot(Iu,sigma_wt(35,:),'-','LineWidth',LW); hold on
plot(Iu,sigma_wt(end,:),'-','LineWidth',LW); hold on

xlabel("$Iu [-]$",'Interpreter','latex','FontSize',FS);
ylabel("$\sigma_{w,t}^2 [m^2]$",'Interpreter','latex','FontSize',FS);
strLegend = {string(round(Lu(1))),string(round(Lu(10))),string(round(Lu(25))),string(round(Lu(35))),string(round(Lu(end)))};
legend("$Lu = " + strLegend + " $",'Interpreter','latex','Location','northwest','FontSize',FS);
grid on;

figure()
plot(Lu,sigma_mbr(:,1)/(10^6),'-','LineWidth',LW); hold on
plot(Lu,sigma_mbr(:,2)/(10^6),'-','LineWidth',LW); hold on
plot(Lu,sigma_mbr(:,3)/(10^6),'-','LineWidth',LW); hold on
plot(Lu,sigma_mbr(:,4)/(10^6),'-','LineWidth',LW); hold on
plot(Lu,sigma_mbr(:,5)/(10^6),'-','LineWidth',LW); hold on

xlabel("$Lu [m]$",'Interpreter','latex','FontSize',FS);
ylabel("$\sigma_{m_b,r}^2 [MN^2 m^2]$",'Interpreter','latex','FontSize',FS);
legend("$Iu = " + string(Iu) + " $",'Interpreter','latex','Location','northwest','FontSize',FS);
grid on;

figure()
plot(Iu,sigma_mbr(1,:)/(10^6),'-','LineWidth',LW); hold on
plot(Iu,sigma_mbr(10,:)/(10^6),'-','LineWidth',LW); hold on
plot(Iu,sigma_mbr(25,:)/(10^6),'-','LineWidth',LW); hold on
plot(Iu,sigma_mbr(35,:)/(10^6),'-','LineWidth',LW); hold on
plot(Iu,sigma_mbr(end,:)/(10^6),'-','LineWidth',LW); hold on

xlabel("$Iu [-]$",'Interpreter','latex','FontSize',FS);
ylabel("$\sigma_{m_b,r}^2 [MN^2 m^2]$",'Interpreter','latex','FontSize',FS);
strLegend = {string(round(Lu(1))),string(round(Lu(10))),string(round(Lu(25))),string(round(Lu(35))),string(round(Lu(end)))};
legend("$Lu = " + strLegend + " $",'Interpreter','latex','Location','northwest','FontSize',FS);
grid on;

%% 
function [Ay,Az,Ayz] = pAxes2aAxes(A1,A2,alpha)
    % According to Bauchau page 234 equations (6.34a,b,c) where for the
    % definition of principal axes we have set A12 = 0
    
    Ay       = 0.5*(A1 + A2) + 0.5*(A1 - A2).*cos(2*alpha);
    Az       = 0.5*(A1 + A2) - 0.5*(A1 - A2).*cos(2*alpha);
    Ayz      = 0.5*(A1 - A2).*sin(2.*alpha);
end