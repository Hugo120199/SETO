clc
clear
close all

load X_2013.mat
load y_2013.mat
load Z_2013.mat

% AUTORREGRESSIVO
% p = 3;
% hp = 1;
% 
% mod_AR = init_AR(p,hp);
% 
% mod_AR_trained = train_AR(mod_AR,y_2013,X_2013,Z_2013);
% 
% y_pred_AR = pred_AR(mod_AR_trained,y_2013,X_2013,Z_2013);

% RED NEURONAL
p = 20;
hp = 12;
estruc = [16,8,4];

mod_ANN = init_ANN(p,estruc,hp);

mod_ANN_trained = train_ANN(mod_ANN,y_2013,X_2013,Z_2013);

y_pred_ANN = pred_ANN(mod_ANN_trained,y_2013,X_2013,Z_2013);


% PLOTS
figure(1)
plot(y_2013); hold on;
%plot(y_pred_AR); hold on;
plot(y_pred_ANN); hold on;
legend('Datos','Predicción')

grid on

% Errores

%Err_AR = eval_RMSE(y_2013,y_pred_AR)
Err_ANN = eval_RMSE(y_2013,y_pred_ANN)


