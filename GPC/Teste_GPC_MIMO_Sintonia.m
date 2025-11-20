%% ========================================================================
% Teste_GPC_MIMO_Sintonia.m
% Teste em malha fechada do GPC da coluna de destilação
% Usando modelo linear MIMO (Gp_d) + perturbação medida (Gd_d)
% ========================================================================
run('GPC_MIMO_Prediction_Matrices.m');
clear; close all; clc;

%% Carregar controlador
% Certifique-se de ter rodado 'GPC_MIMO_Prediction_Matrices.m' (o novo)
load('GPC_MIMO_Controller_Data.mat','controller_data');

Ts  = controller_data.Ts;
N2  = controller_data.N2;
Nu  = controller_data.Nu;

% ESTE É O MODELO POSICIONAL (PLANTA REAL)
A   = controller_data.A;
B   = controller_data.B;
C   = controller_data.C;
B_v = controller_data.B_v;

ny  = controller_data.ny;
nu  = controller_data.nu;
nv  = controller_data.nv;
nx  = size(A,1); % nx da planta real (posicional)

u0  = controller_data.u0;
y0  = controller_data.y0;
v0  = controller_data.v0;

%% Tempo de simulação
T_sim = 250;          % número de amostras
t = (0:T_sim-1)' * Ts;

%% Alocação de variáveis
y = zeros(ny, T_sim);
u = zeros(nu, T_sim);
v = zeros(nv, T_sim);
r = zeros(ny, T_sim);
x = zeros(nx, T_sim);

%% Condições iniciais
y(:,1) = y0;
u(:,1) = u0;
v(:,1) = v0;
x(:,1) = zeros(nx,1);      % estado de referência (desvio zero)

% Referências: degrau em xD e xB
r(:,1:20)   = repmat(y0, 1, 20);                    % regime
r(:,21:end) = repmat([0.95; 0.03], 1, T_sim-20);    % pequenos degraus

% Perturbação medida: degrau em F a partir de k=100
for k = 1:T_sim
    if k < 100
        v(:,k) = v0;
    else
        v(:,k) = v0 + [0.1*sin(T_sim); -0.1*sin(T_sim)];   % exemplo F e qF mudam um pouco
    end
end

%% Loop de simulação
fprintf('Iniciando simulação do GPC...\n');

exitflags = zeros(1,T_sim);

for k = 2:T_sim
    
    % Chamar controlador (usa y(k-1), u(k-1), v(k-1), r(k))
    [u_k, fval, exitflag] = GPC_MIMO_Controller( ...
        y(:,k-1), u(:,k-1), v(:,k-1), r(:,k), controller_data);
    
    u(:,k) = u_k;
    exitflags(k) = exitflag;
    
    % Simular planta linear em desvio do ponto de operação
    u_t = u(:,k) - u0;
    v_t = v(:,k) - v0;
    
    x(:,k) = A*x(:,k-1) + B*u_t  + B_v*v_t;
    y_t    = C*x(:,k);
    y(:,k) = y0 + y_t;
end

fprintf('Simulação concluída. exitflags únicos: '); 
disp(unique(exitflags));

%% Plots

figure('Position', [100, 100, 1280, 720], 'Color', 'White');
subplot(3,1,1);
plot(t, y(1,:), 'b', t, r(1,:), 'r--','LineWidth',1.5);
ylim([0.9, 1]);
grid on; ylabel('x_D');
legend('x_D','r_D', 'Location', 'best', 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
title("Resposta do sistema ao controlador GPC");
ax = gca;
ax.Color = 'white';
ax.XColor = 'black';
ax.YColor = 'black';
ax.GridColor = [0.3 0.3 0.3];
ax.GridAlpha = 0.7;

subplot(3,1,2);
plot(t, y(2,:), 'b', t, r(2,:), 'r--','LineWidth',1.5);
ylim([0.0, 0.1]);
grid on; ylabel('x_B');
legend('x_B','r_B', 'Location', 'best', 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
ax = gca;
ax.Color = 'white';
ax.XColor = 'black';
ax.YColor = 'black';
ax.GridColor = [0.3 0.3 0.3];
ax.GridAlpha = 0.7;

subplot(3,1,3);
stairs(t, u(1,:), 'b','LineWidth',1.5); hold on;
stairs(t, u(2,:), 'r','LineWidth',1.5);
grid on; ylabel('u'); xlabel('Tempo (min)');
legend('V','R', 'Location', 'best', 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
ax = gca;
ax.Color = 'white';
ax.XColor = 'black';
ax.YColor = 'black';
ax.GridColor = [0.3 0.3 0.3];
ax.GridAlpha = 0.7;

% Aplicar configurações de fonte para todos os eixos
allAxes = findobj(gcf, 'Type', 'axes');
for i = 1:length(allAxes)
    set(allAxes(i), 'FontSize', 12, 'FontName', 'Arial');
    
    % Garantir que títulos e labels tenham cor preta
    titleHandle = get(allAxes(i), 'Title');
    xlabelHandle = get(allAxes(i), 'XLabel');
    ylabelHandle = get(allAxes(i), 'YLabel');
    set(titleHandle, 'Color', 'black');
    set(xlabelHandle, 'Color', 'black');
    set(ylabelHandle, 'Color', 'black');
end