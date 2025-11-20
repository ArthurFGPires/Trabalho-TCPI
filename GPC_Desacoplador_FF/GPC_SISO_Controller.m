function [u, fval, exitflag] = GPC_SISO_Controller(y, u_prev, r, controller_data)
% Controlador GPC Monovariável (SISO)

%% Extrair parâmetros
A   = controller_data.A;
B   = controller_data.B;
C   = controller_data.C;
% B_v é ignorado
nx  = controller_data.nx; % nx_aug
nx_pos = controller_data.nx_pos; % nx original

G   = controller_data.G;
F   = controller_data.F;
% H é ignorado

Hess      = controller_data.Hess;
grad_base = controller_data.grad_base;
constraints = controller_data.constraints;

N2  = controller_data.N2;
Nu  = controller_data.Nu;
ny  = controller_data.ny; % ny = 1
nu  = controller_data.nu; % nu = 1

u0  = controller_data.u0;
y0  = controller_data.y0;

%% Inicialização persistentes
persistent x_a_prev u_prev_prev DU_prev
% x_a_prev = estado aumentado x_a(k-1)
% u_prev_prev = u(k-2)
% DU_prev = \Delta U ótimo da iteração anterior

if isempty(x_a_prev)
    x_a_prev = zeros(nx, 1);
end
if isempty(u_prev_prev)
    u_prev_prev = u0;
end
if isempty(DU_prev)
    DU_prev = zeros(nu*Nu, 1);
end

%% Coordenadas de desvio (agora tudo escalar)
y_t_med = y - y0;
u_prev_t = u_prev - u0;
u_prev_prev_t = u_prev_prev - u0;
r_t = r - y0;

%% 1. Calcular \Delta u(k-1)
du_prev_t = u_prev_t - u_prev_prev_t;

%% 2. Observador de Estado e Distúrbio
y_t_pred = x_a_prev(nx_pos+1 : end);
dist_k_minus_1 = y_t_med - y_t_pred;
x_a_prev(nx_pos+1 : end) = y_t_med;

%% 3. Propagação do Estado
% x_a(k) = A_a * x_a(k-1) + B_a * \Delta u(k-1)
x_a_k = A * x_a_prev + B * du_prev_t;

%% 4. Vetores futuros
r_future = repmat(r_t, N2, 1);
v_future_delta = zeros(0, 1); % Sem feedforward

%% 5. Predição livre
y_free_raw = F * x_a_k; % H removido
y_free_t = y_free_raw + repmat(dist_k_minus_1, N2, 1);

%% 6. Gradiente
grad = grad_base * (r_future - y_free_t);

%% 7. Restrições
C_du = constraints.C_du;
d_du = constraints.d_du;

C_u   = constraints.C_u;
Lmat  = constraints.L;
u_min = constraints.u_min;
u_max = constraints.u_max;

u_max_rep = repmat(u_max - u0, Nu, 1);
u_min_rep = repmat(u_min - u0, Nu, 1);

d_u_upper = u_max_rep - Lmat * u_prev_t;
d_u_lower = -(u_min_rep - Lmat * u_prev_t);
d_u = [d_u_upper; d_u_lower];

% Removendo restrições de saída para garantir factibilidade
C_ineq = [C_du; C_u];
d_ineq = [d_du; d_u];

%% 8. Warm start
x0 = DU_prev; 

%% 9. QP
options = optimoptions('quadprog', ...
    'Algorithm','active-set', ...
    'Display','off');

[DU_opt, fval, exitflag] = quadprog( ...
    Hess, grad, ...
    C_ineq, d_ineq, ...
    [], [], [], [], ...
    x0, options);

%% Diagnóstico se falha
if exitflag ~= 1
    warning('GPC-SISO: QP não convergiu - exitflag=%d', exitflag);
    DU_opt = zeros(nu*Nu,1);
else
    DU_prev = DU_opt;
end

%% 10. Primeira ação de controle
Du_t = DU_opt(1); % \Delta u_t(k) (agora é escalar)
u_t  = u_prev_t + Du_t;
u    = u0 + u_t;        

u = max(u_min + u0, min(u_max + u0, u));

%% 11. Atualizar persistentes
x_a_prev = x_a_k;
u_prev_prev = u_prev; 

end