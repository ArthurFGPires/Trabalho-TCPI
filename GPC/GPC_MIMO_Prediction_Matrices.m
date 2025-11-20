%% ========================================================================
% GPC_MIMO_Prediction_Matrices.m
% Cálculo das matrizes para GPC MIMO (CARIMA) com Restrições
% ========================================================================

clear; clc; close all;

%% 1. Definição do Modelo Contínuo
Ts = 1; % min
s = tf('s');

% --- Modelo do Processo Gp(s) (V, R -> xD, xB) ---
G11 = -1.2/(75.6*s + 1) * exp(-0.16*s);   
G12 = 0.4/(54.2*s + 1) * exp(-0.156*s);   
G21 = -0.5/(58*s + 1) * exp(-0.16*s);     
G22 = 1.3/(76.05*s + 1);                  

Gp_c = [G11 G12; G21 G22];

% --- Modelo de Perturbação Gd(s) (F, qF -> xD, xB) ---
Gd_F_xD = 0.12/(52.2*s + 1);
Gd_F_xB = 0.7 /(35.26*s + 1);
Gd_qF_xD = 0.77/(73.2*s + 1) * exp(-0.285*s);
Gd_qF_xB = 1.1 /(72.58*s + 1);

Gd_c = [ Gd_F_xD   Gd_qF_xD;
         Gd_F_xB   Gd_qF_xB ];

%% 2. Discretização

% 2.1 Discretizar Planta
sys_p_d = c2d(Gp_c, Ts, 'zoh');
[Ap, Bp, Cp, ~] = ssdata(sys_p_d);

% 2.2 Discretizar Perturbação
sys_d_d = c2d(Gd_c, Ts, 'zoh');
[Ad, Bd, Cd, ~] = ssdata(sys_d_d);

% Dimensões auxiliares
nx_p = size(Ap, 1);
nx_d = size(Ad, 1);
ny   = size(Cp, 1); % deve ser igual a size(Cd, 1)
nu   = size(Bp, 2);
nv   = size(Bd, 2);

% 2.3 Fusão dos Modelos (Concatenação)
% Sistema "diagonal" onde os estados não se misturam fisicamente,
% mas o controlador enxerga tudo como um único sistema linear.

% A_total: [Ap  0 ]
%          [0   Ad]
A = blkdiag(Ap, Ad);

% B_total: [Bp] (Afeta apenas estados da planta)
%          [0 ]
B = [Bp; zeros(nx_d, nu)];

% B_v_total: [0 ] (Afeta apenas estados da perturbação)
%            [Bd]
B_v = [zeros(nx_p, nv); Bd];

% C_total: [Cp  Cd] (A saída é a soma dos efeitos: y = y_p + y_d)
C = [Cp, Cd];

nx = size(A, 1); % nx total = nx_p + nx_d

%% 3. Parâmetros do GPC
N2 = 20;                % Horizonte de Predição
Nu = 15;                % Horizonte de Controle
Lambda = 5e4;           % Peso para ação de controle
Qy = diag([100, 100]);  % Peso para erro de referência
Ru = Lambda * eye(nu);  % Matriz diagonal com peso da ação de controle

%% 4. Cálculo das Matrizes Dinâmicas (G e H)
% O cálculo via resposta ao degrau

G = zeros(ny*N2, nu*Nu);
H = zeros(ny*N2, nv*N2);

% --- Construção de G (Step Response de u) ---
StepCoeffs_u = zeros(ny, nu, N2);
for j = 1:nu
    x_sim = zeros(nx, 1);
    u_step = zeros(nu, 1); u_step(j) = 1;
    for k = 1:N2
        x_sim = A*x_sim + B*u_step;
        StepCoeffs_u(:, j, k) = C*x_sim;
    end
end

for i = 1:N2
    for j = 1:min(i, Nu)
        G( (i-1)*ny+1 : i*ny, (j-1)*nu+1 : j*nu ) = StepCoeffs_u(:, :, i-j+1);
    end
end

% --- Construção de H (Step Response de v) ---
StepCoeffs_v = zeros(ny, nv, N2);
for j = 1:nv
    x_sim = zeros(nx, 1);
    v_step = zeros(nv, 1); v_step(j) = 1;
    for k = 1:N2
        x_sim = A*x_sim + B_v*v_step;
        StepCoeffs_v(:, j, k) = C*x_sim;
    end
end
for i = 1:N2
    for j = 1:i
        H( (i-1)*ny+1 : i*ny, (j-1)*nv+1 : j*nv ) = StepCoeffs_v(:, :, i-j+1);
    end
end

%% 5. Matrizes do QP
Q_bar = kron(eye(N2), Qy);
R_bar = kron(eye(Nu), Ru);

Hess = 2 * (G' * Q_bar * G + R_bar);            % Definindo matriz Hessiana
Hess = (Hess + Hess')/2 + 1e-6*eye(size(Hess)); % Normalizando (evitar flag do QP)
grad_base = -2 * G' * Q_bar;                    % Cálculo do gradiente do termo linear

%% 6. Restrições (Entrada e Saída)

u_abs_min = [0; 0];     % Ação de controle minima
u_abs_max = [12; 10];   % Ação de controle máxima
du_lim    = [1; 1]; % Máxima variação da ação de controle entre passos

u0 = [2.4; 1.9];        % Ação de controle dos pontos de operação
y0 = [0.96; 0.04];      % Saídas do ponto de operação
v0 = [1.0; 1.0];        % Valor das perutrbações do ponto de operação

% Matrizes para Delta U
C_du = [ eye(nu*Nu); -eye(nu*Nu) ];                     
d_du = [ repmat(du_lim, Nu, 1); repmat(du_lim, Nu, 1) ];

% Matrizes para U
E = zeros(nu*Nu, nu*Nu);

for i = 1:Nu
    for j = 1:i
        E( (i-1)*nu+1:i*nu, (j-1)*nu+1:j*nu ) = eye(nu);
    end
end

C_u = [ E; -E ];
L = kron(ones(Nu,1), eye(nu));

% Restrições de Saída
y_abs_min = [0; 0];    % valor mínimo permitido as saídas 
y_abs_max = [1; 1];    % valor máximo permitido as saídas
C_y = [ G; -G ];

% Salvar dados das restrições
constraints.C_du  = C_du;
constraints.d_du  = d_du;

constraints.C_u   = C_u;
constraints.L     = L;
constraints.u_abs_min = u_abs_min;
constraints.u_abs_max = u_abs_max;

constraints.C_y   = C_y;
constraints.y_abs_min = y_abs_min;
constraints.y_abs_max = y_abs_max;

%% 7. Salvar Dados do controlador
controller_data.A = A; controller_data.B = B; controller_data.C = C; controller_data.B_v = B_v;
controller_data.nx = nx;
controller_data.G = G; controller_data.H = H;
controller_data.Hess = Hess; controller_data.grad_base = grad_base;
controller_data.constraints = constraints;
controller_data.N2 = N2; controller_data.Nu = Nu;
controller_data.ny = ny; controller_data.nu = nu; controller_data.nv = nv;
controller_data.Ts = Ts;
controller_data.u0 = u0; controller_data.y0 = y0; controller_data.v0 = v0;

save('GPC_MIMO_Controller_Data.mat', 'controller_data');
fprintf('Dados salvos.\n');