%% ========================================================================
% RUN_Build_Multiloop_GPCs.m
% Constrói os dois controladores SISO para a arquitetura Multiloop
% Também define os parâmetros da planta 'p' para a simulação.
% ========================================================================

clear; clc; close all;

%% ----------------- Definição dos Parâmetros da Planta (p) -----------------
% Parâmetros físicos da coluna de destilação, conforme 'Trabalho_TCPI.pdf'

p.NT   = 28;     % Número de estágios (pratos) [cite: 56]
p.NF   = 15;     % Estágio de alimentação [cite: 57]
p.alpha    = 1.6;    % Volatilidade relativa (alpha) [cite: 58]
p.M0   = 1.0;    % Retenção de líquido nominal (M0) [kmol] [cite: 59]
p.tau_h   = 2.0;    % Constante de tempo hidráulica (tau_h) [min] [cite: 60]

% Ponto de Operação Nominal (usado na simulação)
p.F   = 1.0;    % Vazão de alimentação [Kmol/min] [cite: 61]
p.zF  = 0.5;    % Composição da alimentação [cite: 62]
p.qF  = 1.0;    % Qualidade da alimentação (líquido saturado) [cite: 63]
p.xD  = 0.96;   % Composição de topo desejada [cite: 65]
p.xB  = 0.04;   % Composição de fundo desejada [cite: 66]
p.V   = 2.4;    % Vazão de vapor nominal [kmol/min] [cite: 168]
p.R   = 1.9;    % Vazão de refluxo nominal [kmol/min] [cite: 168]
x0ss = linspace(0.01, 0.99, p.NT);
M0_ss = p.M0;

%% ----------------- Definições da Planta (Modelos) -----------------
Ts = 10;
s = tf('s');

% Modelos (conforme GPC_MIMO_Prediction_Matrices.m)
% Usando o emparelhamento do TCC (V->xD, R->xB)
G11 = -1.2/(75.6*s + 1) * exp(-0.16*s);   % V -> xD [cite: 318]
G22 = 1.3/(76.05*s + 1);                  % R -> xB [cite: 307]

% Pontos de operação (já definidos em 'p')
V0 = p.V;  xD0 = p.xD;
R0 = p.R;  xB0 = p.xB;

%% ----------------- Sintonia (Pode ajustar independentemente) -----------------
% Sintonia Loop 1 (V -> xD)
N2_1 = 40; Nu_1 = 15;
Qy_1 = 100; Lambda_1 = 50; % Ex: Lambda alto para acalmar Loop 1

% Sintonia Loop 2 (R -> xB)
N2_2 = 40; Nu_2 = 15;
Qy_2 = 100; Lambda_2 = 50; % Ex: Lambda alto para acalmar Loop 2

%% ----------------- Gerar Controladores -----------------
% Gerar Loop 1
GPC_SISO_Matrices(G11, Ts, V0, xD0, N2_1, Nu_1, Qy_1, Lambda_1, 'GPC_Loop1_Data.mat');

% Gerar Loop 2
GPC_SISO_Matrices(G22, Ts, R0, xB0, N2_2, Nu_2, Qy_2, Lambda_2, 'GPC_Loop2_Data.mat');

fprintf('=== Controladores Multiloop gerados! ===\n');
fprintf('=== Estrutura de parâmetros ''p'' criada no workspace. ===\n');

% O próximo passo seria chamar seu script de simulação:
model_name = 'coluna_sim_multiloop';
disp('Abrindo modelo Simulink e iniciando simulação...');
open_system(model_name);
sim(model_name);
disp('Simulação concluída.');
