%% Modelos de predicción de potencia
clc 
clear all
close all 

%% Carga de datos de las series temporales de entrenamiento
load X_2013.mat
load y_2013.mat
load Z_2013.mat


%% Modelo Autorregresivo (AR) 
hp = 1;

    % Algunos datos del histórico
number = 1000:500:8000;
number = [number length(y_2013)]; 
number = 2145:5:2165; 

for j = 1:length(number)-1 
    y_2013 = y_2013(1,end-number(j):end);
    X_2013 = X_2013(1,end-number(j):end);
    Z_2013 = Z_2013(2,end-number(j):end);
    
    p = 1:1:100;
    for i = 1:length(p)
        modHp1_AR = init_AR(p(i),hp);
        modHp1_AR_trained = train_AR(modHp1_AR,y_2013,0,0);
        y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,0,0);
        ErrY_AR_hp1(i,j) = eval_RMSE(y_2013,y_predHp1_AR);
    
        modHp1_AR = init_AR(p(i),hp);
        modHp1_AR_trained = train_AR(modHp1_AR,y_2013,X_2013,0);
        y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,X_2013,0);
        ErrYX_AR_hp1(i,j) = eval_RMSE(y_2013,y_predHp1_AR);
    
        modHp1_AR = init_AR(p(i),hp);
        modHp1_AR_trained = train_AR(modHp1_AR,y_2013,0,Z_2013);
        y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,0,Z_2013);
        ErrYZ_AR_hp1(i,j) = eval_RMSE(y_2013,y_predHp1_AR);
    
        modHp1_AR = init_AR(p(i),hp);
        modHp1_AR_trained = train_AR(modHp1_AR,y_2013,X_2013,Z_2013);
        y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,X_2013,Z_2013);
        ErrYXZ_AR_hp1(i,j) = eval_RMSE(y_2013,y_predHp1_AR);
    end
    
    figure(j)
    plot(p,ErrY_AR_hp1(:,j),'-r'); hold on; 
    plot(p,ErrYX_AR_hp1(:,j),'-b'); hold on; 
    plot(p,ErrYZ_AR_hp1(:,j),'-g'); hold on; 
    plot(p,ErrYXZ_AR_hp1(:,j),'-m'); hold on; 
    grid on;
    xlabel('p','interpreter','latex','fontsize',16)
    ylabel('RMSE','interpreter','latex','fontsize',16)
    legend('Hist\''orico Y', 'Hist\''orico Y,X','Hist\''orico Y,Z',...
        'Hist\''orico Y,X,Z','interpreter','latex','fontsize',16)
    sgtitle([num2str(number(j)),' valores del hist\''orico'],...
        'interpreter','latex','fontsize',16)
    
    load X_2013.mat
    load y_2013.mat
    load Z_2013.mat
end 

    % Todos los datos del histórico
p = 1000:10:2000;
j = 0;
for i = 1:length(p)
    modHp1_AR = init_AR(p(i),hp);
    modHp1_AR_trained = train_AR(modHp1_AR,y_2013,0,0);
    y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,0,0);
    ErrY_AR_hp1(i,j+1) = eval_RMSE(y_2013,y_predHp1_AR);

    modHp1_AR = init_AR(p(i),hp);
    modHp1_AR_trained = train_AR(modHp1_AR,y_2013,X_2013,0);
    y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,X_2013,0);
    ErrYX_AR_hp1(i,j+1) = eval_RMSE(y_2013,y_predHp1_AR);

    modHp1_AR = init_AR(p(i),hp);
    modHp1_AR_trained = train_AR(modHp1_AR,y_2013,0,Z_2013);
    y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,0,Z_2013);
    ErrYZ_AR_hp1(i,j+1) = eval_RMSE(y_2013,y_predHp1_AR);

    modHp1_AR = init_AR(p(i),hp);
    modHp1_AR_trained = train_AR(modHp1_AR,y_2013,X_2013,Z_2013);
    y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,X_2013,Z_2013);
    ErrYXZ_AR_hp1(i,j+1) = eval_RMSE(y_2013,y_predHp1_AR);
end
figure(j+1)
plot(p,ErrY_AR_hp1(:,j+1),'-r'); hold on; 
plot(p,ErrYX_AR_hp1(:,j+1),'-b'); hold on; 
plot(p,ErrYZ_AR_hp1(:,j+1),'-g'); hold on; 
plot(p,ErrYXZ_AR_hp1(:,j+1),'-m'); hold on; 
grid on;
xlabel('Ventana temporal, p','interpreter','latex','fontsize',16)
ylabel('Error medio cuadr\''atico, RMSE','interpreter','latex','fontsize',16)
legend('Hist\''orico Y', 'Hist\''orico Y,X','Hist\''orico Y,Z',...
    'Hist\''orico Y,X,Z','interpreter','latex','fontsize',16)
sgtitle([num2str(length(y_2013)),' valores del hist\''orico'],...
    'interpreter','latex','fontsize',16)

%% Representación gráfica y cálculo del error
load X_2013.mat
load y_2013.mat
load Z_2013.mat

    % Elegimos una combinación y representamos los resultados de predicción
    % frente a los datos de entrenamiento: 
%       - Nº datos del histórico empleados para entrenar: 2145-2155 [N]
%       - Datos histótico: X, Y y Z (todos en combinación) 
%       - Ventana temporal: 50 [P]

P = 50; 
N = 2150;

y_ = y_2013(1,end-N:end);
X_ = X_2013(1,end-N:end);
Z_ = Z_2013(2,end-N:end);

modHp1_AR = init_AR(P,hp);
modHp1_AR_trained = train_AR(modHp1_AR,y_,X_,Z_);
y_predHp1_AR = pred_AR(modHp1_AR_trained,y_,X_,Z_);
% y_predHp1_AR = pred_AR(modHp1_AR_trained,y_2013,X_2013,Z_2013);
ErrYXZ_AR_hp1 = eval_RMSE(y_,y_predHp1_AR);


    % Representación de las series temporales
figure()
plot(y_2013,'-r'); hold on;
% plot(y_,'-r'); hold on;
plot(y_predHp1_AR,'-b'); hold on;
grid on
xlabel('Medida')
ylabel('Potencia')

