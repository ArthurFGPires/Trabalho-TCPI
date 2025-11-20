function [u, fval, exitflag] = GPC_MIMO_Controller(y, u_prev, v, r, controller_data)
% GPC MIMO (CARIMA) com Restrições

%% Extrair parâmetros do controlador
A   = controller_data.A;
B   = controller_data.B;
C   = controller_data.C;
B_v = controller_data.B_v;

Hess      = controller_data.Hess;
grad_base = controller_data.grad_base;
constraints = controller_data.constraints;

N2  = controller_data.N2;
Nu  = controller_data.Nu;
nx  = controller_data.nx;
ny  = controller_data.ny;
nu  = controller_data.nu;

u0  = controller_data.u0;
y0  = controller_data.y0;
v0  = controller_data.v0;

%% Inicialização persistentes
persistent x_est DU_prev v_prev u_prev_mem

if isempty(x_est), x_est = zeros(nx, 1); end
if isempty(DU_prev), DU_prev = zeros(nu*Nu, 1); end
if isempty(v_prev), v_prev = v0; end
if isempty(u_prev_mem), u_prev_mem = u0; end

%% Variáveis de Desvio
y_t = y - y0;
u_prev_t = u_prev - u0;
v_t = v - v0;           
r_t = r - y0;

v_prev_t = v_prev - v0;
u_prev_mem_t = u_prev_mem - u0;

%% 1. Observador de Estados
x_pred = A * x_est + B * u_prev_mem_t + B_v * v_prev_t;
y_model = C * x_pred;
dist_k = y_t - y_model; % Erro de modelagem
x_est = x_pred;

%% 2. Predição Livre (Feedforward Dinâmico)
y_free = zeros(ny * N2, 1);
x_sim = x_est;
u_constant = u_prev_t; 
v_constant = v_t;      

for k = 1:N2
    x_sim = A * x_sim + B * u_constant + B_v * v_constant;
    y_free((k-1)*ny + 1 : k*ny) = C * x_sim + dist_k;
end

%% 3. Otimização QP
r_future = repmat(r_t, N2, 1);
grad = grad_base * (r_future - y_free);

% --- Restrições de Delta U (Taxa) ---
C_du = constraints.C_du;
d_du = constraints.d_du;

% --- Restrições de U (Amplitude) ---
C_u   = constraints.C_u;
Lmat  = constraints.L;

% Limites absolutos de entrada convertidos para desvio
u_min_dev = constraints.u_abs_min - u0;
u_max_dev = constraints.u_abs_max - u0;

u_max_rep = repmat(u_max_dev, Nu, 1);
u_min_rep = repmat(u_min_dev, Nu, 1);

% d_u define quanto ainda pode subir ou descer baseado no u_prev atual
d_u_upper = u_max_rep - Lmat * u_prev_t;
d_u_lower = -(u_min_rep - Lmat * u_prev_t);
d_u = [d_u_upper; d_u_lower];

% --- Restrições de Y (Saída) ---
C_y = constraints.C_y;

% Limites absolutos de saída convertidos para desvio
% y_min <= y_total <= y_max  -->  y_min - y0 <= y_dev <= y_max - y0
y_min_dev = constraints.y_abs_min - y0;
y_max_dev = constraints.y_abs_max - y0;

y_max_rep = repmat(y_max_dev, N2, 1);
y_min_rep = repmat(y_min_dev, N2, 1);

% A inequação é: G*du <= y_max_dev - y_free
%               -G*du <= -(y_min_dev - y_free)
d_y_upper = y_max_rep - y_free;
d_y_lower = -(y_min_rep - y_free);
d_y = [d_y_upper; d_y_lower];

% --- Montagem Final ---
C_ineq = [C_du; C_u; C_y];
d_ineq = [d_du; d_u; d_y];

% Solver
options = optimoptions('quadprog', 'Algorithm','active-set', 'Display','off');
[DU_opt, fval, exitflag] = quadprog(Hess, grad, C_ineq, d_ineq, [], [], [], [], DU_prev, options);

if exitflag ~= 1
    DU_opt = zeros(nu*Nu, 1);
    fprintf('QP Falhou! Exitflag: %f\n', exitflag);
end
DU_prev = DU_opt;

%% 4. Aplicação do Controle
du_k = DU_opt(1:nu);
u_t = u_prev_t + du_k;
u = u_t + u0;

% Saturação final de segurança
u = max(constraints.u_abs_min, min(constraints.u_abs_max, u));

%% 5. Atualização de Memória
u_prev_mem = u;
v_prev = v; 

end