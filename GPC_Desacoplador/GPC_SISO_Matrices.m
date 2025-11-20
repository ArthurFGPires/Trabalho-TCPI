function GPC_SISO_Matrices(Gp_c, Gd_c, Ts, u0_siso, y0_siso, v0_full, N2, Nu, Qy, Lambda, save_filename)
% GPC_SISO_Matrices.m
% Calcula matrizes para GPC SISO com Feedforward Dinâmico (Abordagem CARIMA)
%
% Gp_c: TF contínua do processo (1x1)
% Gd_c: TF contínua das perturbações (1xNv)
% v0_full: vetor com os valores nominais de todas as perturbações

%% ----------------- Modelo e Discretização -----------------
% Combina planta e perturbação para discretização coerente de atrasos
sys_c = [Gp_c, Gd_c]; 

% Discretização (ZOH)
sys_d = c2d(sys_c, Ts, 'zoh');

% Extrair as matrizes de estado no formato Posicional
[A, B_total, C, ~] = ssdata(sys_d);

nu = size(Gp_c, 2); % nu = 1 (SISO)
nv = size(Gd_c, 2); % nv (ex: 2 para F e qF)
ny = 1;
nx = size(A, 1);

% Separar a matriz de entrada B_total em B (controle) e B_v (perturbação)
B   = B_total(:, 1:nu);       % Primeira coluna
B_v = B_total(:, nu+1:end);   % Restante das colunas

fprintf('Calculando GPC SISO para %s...\n', save_filename);
fprintf('Modelo Discreto: nx=%d, nu=%d, ny=%d, nv=%d\n', nx, nu, ny, nv);

%% ----------------- Horizontes e pesos -----------------
Ru = Lambda * eye(nu);

%% ----------------- Cálculo de G e H via Resposta ao Degrau ----------------
% G: Efeito de u na saída
% H: Efeito de v na saída (Feedforward)

G = zeros(ny*N2, nu*Nu);
H = zeros(ny*N2, nv*N2);

% --- Construção de G (Step Response de u) ---
% Como nu=1, só temos um canal de controle
StepCoeffs_u = zeros(ny, N2);
x_sim = zeros(nx, 1);
u_step = 1; 

for k = 1:N2
    x_sim = A*x_sim + B*u_step;
    StepCoeffs_u(k) = C*x_sim;
end

for i = 1:N2
    for j = 1:min(i, Nu)
        G(i, j) = StepCoeffs_u(i - j + 1);
    end
end

% --- Construção de H (Step Response de v) ---
% Para cada perturbação (coluna de B_v), calculamos a resposta
StepCoeffs_v = zeros(ny, nv, N2);

for j = 1:nv
    x_sim = zeros(nx, 1);
    v_step = zeros(nv, 1); 
    v_step(j) = 1; % Degrau unitário na perturbação j
    
    for k = 1:N2
        x_sim = A*x_sim + B_v*v_step;
        StepCoeffs_v(1, j, k) = C*x_sim;
    end
end

% Monta H (Blocos triangulares)
for i = 1:N2
    for j = 1:i
        % H espera blocos de tamanho (ny x nv). Aqui ny=1.
        H(i, (j-1)*nv+1 : j*nv) = StepCoeffs_v(1, :, i-j+1);
    end
end

%% ----------------- Matrizes da função custo (QP) ----------------------
Q_bar  = kron(eye(N2), Qy);
R_bar  = kron(eye(Nu), Ru);

Hess   = 2*(G' * Q_bar * G + R_bar);
Hess   = (Hess + Hess') / 2 + 1e-6 * eye(size(Hess)); % Regularização
grad_base = -2*G' * Q_bar;

%% ----------------- Restrições ----------------------
% Definindo limites absolutos e taxas (exemplo genérico, ajuste conforme necessidade)
u_abs_max = 12; 
u_abs_min = 0;
du_max = 12; 
du_min = -12;

% Restrições em Delta U
C_du = [ eye(nu*Nu); -eye(nu*Nu) ];
d_du = [ repmat(du_max, Nu, 1); repmat(-du_min, Nu, 1) ];

% Restrições em U (Matriz triangular de soma)
E = tril(ones(Nu));
C_u = [ E; -E ];
L = ones(Nu,1);

constraints.C_du  = C_du;
constraints.d_du  = d_du;
constraints.C_u   = C_u;
constraints.L     = L;
constraints.u_abs_min = u_abs_min;
constraints.u_abs_max = u_abs_max;

%% ----------------- Salvar dados -----------------
controller_data.A = A;
controller_data.B = B;
controller_data.C = C;
controller_data.B_v = B_v; 
controller_data.nx = nx;

% Mantém campos para compatibilidade
controller_data.A_pos = A;
controller_data.B_pos = B;
controller_data.C_pos = C;
controller_data.nx_pos = nx;

controller_data.G = G; 
controller_data.H = H; % Salvo para referência
controller_data.Hess = Hess; 
controller_data.grad_base = grad_base;
controller_data.constraints = constraints;

controller_data.N2 = N2; 
controller_data.Nu = Nu;
controller_data.ny = ny; 
controller_data.nu = nu; 
controller_data.nv = nv;

controller_data.u0 = u0_siso;
controller_data.y0 = y0_siso;
controller_data.v0 = v0_full; 

save(save_filename,'controller_data');
fprintf('Arquivo salvo: %s (SISO com Feedforward Dinâmico)\n\n', save_filename);
end