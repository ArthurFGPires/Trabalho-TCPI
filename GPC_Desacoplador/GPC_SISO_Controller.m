function [u, fval, exitflag] = GPC_SISO_Controller(y, u_prev, v, r, controller_data)
% Controlador GPC SISO (MISO) - Versão CARIMA com Feedforward Dinâmico
% Sem lógica 'in_manual'.

%% Extrair parâmetros
A   = controller_data.A;
B   = controller_data.B;
C   = controller_data.C;
B_v = controller_data.B_v;
nx  = controller_data.nx; 

Hess      = controller_data.Hess;
grad_base = controller_data.grad_base;
constraints = controller_data.constraints;

N2  = controller_data.N2;
Nu  = controller_data.Nu;
ny  = controller_data.ny;
nu  = controller_data.nu;
nv  = controller_data.nv;

u0  = controller_data.u0;
y0  = controller_data.y0;
v0  = controller_data.v0;

%% Inicialização persistentes
persistent x_est DU_prev v_prev u_prev_mem

if isempty(x_est)
    x_est = zeros(nx, 1); % Estado estimado (desvio)
end
if isempty(DU_prev)
    DU_prev = zeros(nu*Nu, 1);
end
if isempty(v_prev)
    v_prev = v0;
end
if isempty(u_prev_mem)
    u_prev_mem = u0;
end

%% Coordenadas de desvio
y_t = y - y0;
u_prev_t = u_prev - u0;
v_t = v - v0;          
r_t = r - y0;

% Valores do passo anterior (para reconstrução do estado)
v_prev_t = v_prev - v0;
u_prev_mem_t = u_prev_mem - u0;

%% 1. Atualização do Estado (Observador)
% Estima o estado atual (k) com base na evolução do passo anterior (k-1 -> k)
% x(k) = A*x(k-1) + B*u(k-1) + Bv*v(k-1)
x_pred = A * x_est + B * u_prev_mem_t + B_v * v_prev_t;

% Saída predita pelo modelo
y_model = C * x_pred;

% Cálculo do Distúrbio/Erro de Modelagem (Ação Integral)
% d(k) = y_medido(k) - y_modelo(k)
dist_k = y_t - y_model;

% Atualiza estimativa
x_est = x_pred;

%% 2. Predição Livre (Feedforward Dinâmico)
% Simula o modelo para o futuro (k+1 ... k+N2)
% Assume-se Delta u = 0 (entrada constante em u_prev_t)
% O Feedforward ocorre ao usar v_t (valor ATUAL) na propagação.
% A diferença entre v_prev (usado no observador) e v_t (usado aqui) gera a dinâmica.

y_free = zeros(ny * N2, 1);
x_sim = x_est;
u_constant = u_prev_t; 
v_constant = v_t;      

for k = 1:N2
    % Propaga estado
    x_sim = A * x_sim + B * u_constant + B_v * v_constant;
    
    % Saída Livre = Modelo + Distúrbio Estimado (assumido constante)
    y_free(k) = C * x_sim + dist_k;
end

%% 3. Otimização QP
r_future = repmat(r_t, N2, 1);

% Gradiente: grad = grad_base * (r - y_free)
grad = grad_base * (r_future - y_free);

%% 4. Restrições
C_du = constraints.C_du;
d_du = constraints.d_du;
C_u   = constraints.C_u;
Lmat  = constraints.L;

% Limites físicos convertidos para desvio
u_min_local = constraints.u_abs_min - u0;
u_max_local = constraints.u_abs_max - u0;

u_max_rep = repmat(u_max_local, Nu, 1);
u_min_rep = repmat(u_min_local, Nu, 1);

% Restrição de amplitude: u_min <= u_prev + sum(du) <= u_max
d_u_upper = u_max_rep - Lmat * u_prev_t;
d_u_lower = -(u_min_rep - Lmat * u_prev_t);
d_u = [d_u_upper; d_u_lower];

C_ineq = [C_du; C_u];
d_ineq = [d_du; d_u];

%% 5. Solução QP
options = optimoptions('quadprog', 'Algorithm','active-set', 'Display','off');

[DU_opt, fval, exitflag] = quadprog(Hess, grad, C_ineq, d_ineq, [], [], [], [], DU_prev, options);

if exitflag ~= 1
    % Falha: mantém controle anterior (Delta u = 0)
    DU_opt = zeros(nu*Nu,1);
end

DU_prev = DU_opt;

%% 6. Aplicar Controle
Du_t = DU_opt(1);     
u_t  = u_prev_t + Du_t;
u    = u0 + u_t;      % Valor absoluto

% Saturação final de segurança
u = max(constraints.u_abs_min, min(constraints.u_abs_max, u));

%% 7. Atualizar persistentes
u_prev_mem = u; 
v_prev = v; 

end