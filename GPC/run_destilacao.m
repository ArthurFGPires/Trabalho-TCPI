% --- Inicia Dados do controlador GPC -----
run("GPC_MIMO_Prediction_Matrices.m");

% --- Limpa os dados ---
clear; clc;
% --- Parâmetros da Coluna (Cenário 3: Ponto de Operação Específico) ---
% Objetivo: Atingir xD=0.95 e xB=0.05 em malha aberta

p.NT = 28;         % Número total de estágios (aumentado para a separação mais difícil)
p.NF = 15;         % Estágio de alimentação
p.alpha = 1.6;     % Volatilidade relativa (nova especificação)

% --- Condições de Alimentação ---
p.F = 1.0;         % Vazão de alimentação [kmol/min]
p.zF = 0.5;        % Composição da alimentação (50/50 para este projeto)
p.qF = 1.0;        % Qualidade da alimentação (líquido saturado)

% =========================================================================
% --- CONDIÇÕES OPERACIONAIS (PROJETADAS PARA ATINGIR A META) ---
% =========================================================================
% Análise de projeto para alpha=1.6, xD=0.95, xB=0.05:
% - Número Mínimo de Estágios (N_min): ~12.5. Escolhido NT=28.
% - Refluxo Mínimo (R_min): ~2.9.
% - Estratégia: Usar um refluxo operacional R = 1.3 * R_min ≈ 3.8
% - Isso resulta em L = R * D = 3.8 * 0.5 = 1.9
% - E V = L + D = 1.9 + 0.5 = 2.4

p.R = 1.9;         % Vazão de refluxo [kmol/min]
p.V = 2.4;         % Vazão de vapor correspondente [kmol/min]
p.xD = 0.96;
p.xB = 0.04;

% =========================================================================
% --- PARÂMETROS DINÂMICOS (PARA RESPOSTA LENTA E COM ATRASO) ---
% =========================================================================
% Mantém a retenção de líquido grande para uma dinâmica lenta.
p.M0 = 1.0 * ones(p.NT, 1); % Retenção de líquido nominal [kmol]

% Mantém a constante de tempo hidráulica para um atraso relevante.
p.tau_h = 2.0;               % Constante de tempo hidráulica [min]
% =========================================================================

% --- Cálculo do Estado Estacionário Inicial ---
x0_ss = linspace(0.01, 0.99, p.NT)';
M0_ss = p.M0;

% Acelera a chegada ao ponto de operação
x0_ss(1, 1) = 0.04;     
x0_ss(end, 1) = 0.96;   

initial_states = [M0_ss.* x0_ss; M0_ss];

% --- Executar Simulação ---
model_name = 'coluna_sim';
disp('Abrindo modelo Simulink e iniciando simulação...');
open_system(model_name);
sim(model_name);
disp('Simulação concluída.');

% ------ Observações -----------
% Caso receba o erro: Size mismatch for MATLAB expression '<output of load>.controller_data.G'. Expected = 80x30, actual = 80x80.
% É neceesário realizar os seguintes passos
% Delete a pasta 'slprj' e os arquivos 'coluna_sim.slxc', 'GPC_MIMO_Controller_Data.mat'
% Limpe o ambiente com 'clear all'
% Execute o script de cálculo de matrizes e em seguida, sem limpar o ambiente, todos os comandoas até initial_states deste arquivo
% Por ultimo, abra o modelo simulink e pressione Ctrl+D, ao fim da compilação execute sim('coluna_sim.slx')