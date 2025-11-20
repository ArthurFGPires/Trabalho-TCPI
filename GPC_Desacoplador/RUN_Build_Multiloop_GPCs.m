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

% Modelos de Processo (Diagonal)
G11 = -1.2/(75.6*s + 1) * exp(-0.16*s);   % V -> xD
G22 = 1.3/(76.05*s + 1);                  % R -> xB

% *** ADICIONADO: Modelos de Perturbação (do PDF) ***
% Efeito de F (v1) e qF (v2) em xD (y1)
Gd_F_xD = 0.12/(52.2*s + 1);
Gd_qF_xD = 0.77/(73.2*s + 1) * exp(-0.285*s);
Gd_1 = [Gd_F_xD, Gd_qF_xD]; % Linha 1 da matriz Gd

% Efeito de F (v1) e qF (v2) em xB (y2)
Gd_F_xB = 0.7 /(35.26*s + 1);
Gd_qF_xB = 1.1 /(72.58*s + 1);
Gd_2 = [Gd_F_xB, Gd_qF_xB]; % Linha 2 da matriz Gd

% Pontos de operação
V0 = p.V;  xD0 = p.xD;
R0 = p.R;  xB0 = p.xB;
v0 = [p.F; p.qF]; % Vetor de perturbação nominal

%% ----------------- Sintonia (Pode ajustar independentemente) -----------------
% Sintonia Loop 1 (V -> xD)
N2_1 = 40; Nu_1 = 15;
Qy_1 = 100; Lambda_1 = 50; 

% Sintonia Loop 2 (R -> xB)
N2_2 = 40; Nu_2 = 15;
Qy_2 = 100; Lambda_2 = 50; 

%% ----------------- Gerar Controladores -----------------
% Gerar Loop 1 (passando o modelo de perturbação Gd_1)
GPC_SISO_Matrices(G11, Gd_1, Ts, V0, xD0, v0, N2_1, Nu_1, Qy_1, Lambda_1, 'GPC_Loop1_Data.mat');

% Gerar Loop 2 (passando o modelo de perturbação Gd_2)
GPC_SISO_Matrices(G22, Gd_2, Ts, R0, xB0, v0, N2_2, Nu_2, Qy_2, Lambda_2, 'GPC_Loop2_Data.mat');

fprintf('=== Controladores Multiloop (MISO) gerados! ===\n');
fprintf('=== Estrutura de parâmetros ''p'' criada no workspace. ===\n');

% O próximo passo seria chamar seu script de simulação:
model_name = 'coluna_sim_multiloop';
disp('Abrindo modelo Simulink e iniciando simulação...');
open_system(model_name);
sim(model_name);
disp('Simulação concluída.');
