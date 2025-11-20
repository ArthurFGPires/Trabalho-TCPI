function GPC_SISO_Matrices(Gp_c, Ts, u0_siso, y0_siso, N2, Nu, Qy, Lambda, save_filename)
% Calcula as matrizes de predição para UM loop GPC SISO.

%% ----------------- Modelo -----------------
Gp_d = c2d(Gp_c, Ts, 'zoh');
[A_pos, B_pos, C_pos, ~] = ssdata(ss(Gp_d));

ny = 1; nu = 1; nv = 0; % SISO, sem perturbações medidas
nx_pos = size(A_pos, 1);

fprintf('Calculando GPC SISO para %s...\n', save_filename);
fprintf('Modelo POSICIONAL: nx = %d\n', nx_pos);

%% ----------------- Modelo Aumentado (para ação integral) -----------------
A_a = [ A_pos,     zeros(nx_pos, ny);
        C_pos*A_pos, eye(ny)          ];
B_a = [ B_pos;
        C_pos*B_pos ];
C_a = [ zeros(ny, nx_pos),  eye(ny) ];
nx_aug = size(A_a, 1);
fprintf('Modelo AUMENTADO: nx_a = %d\n', nx_aug);

A = A_a; B = B_a; C = C_a; nx = nx_aug;
B_v = zeros(nx_aug, 0); % Sem matriz Bv

%% ----------------- Horizontes e pesos -----------------
Ru = Lambda * eye(nu); % Ru = Lambda (escalar)

%% ----------------- Construção de G, F, H -------------------
G = zeros(ny*N2, nu*Nu);
F = zeros(ny*N2, nx);
H = zeros(ny*N2, nv*N2); % H ficará vazio

for i = 1:N2
    F( (i-1)*ny+1 : i*ny, : ) = C * (A^i);
    % H loop é pulado
    for j = 1:min(i, Nu)
        G( (i-1)*ny+1 : i*ny, (j-1)*nu+1 : j*nu ) = C * (A^(i-j)) * B;
    end
end

%% ----------------- Matrizes da função custo (QP) ----------------------
Q_bar  = kron(eye(N2), Qy);
R_bar  = kron(eye(Nu), Ru);
Hess   = 2*(G' * Q_bar * G + R_bar);
Hess   = (Hess + Hess') / 2 + 1e-6 * eye(size(Hess));
grad_base = -2*G' * Q_bar;

%% ----------------- Restrições ----------------------
% (Valores de exemplo, ajuste se necessário)
u_max = 12 - u0_siso; u_min = 0 - u0_siso;
du_max = 12; du_min = -12;

C_du = [ eye(nu*Nu); -eye(nu*Nu) ];
d_du = [ repmat(du_max, Nu, 1); repmat(-du_min, Nu, 1) ];

E = tril(ones(Nu)); % Simplificado para SISO
C_u = [ E; -E ];
L = ones(Nu,1);

% (Restrições de saída removidas)
constraints.C_du  = C_du;
constraints.d_du  = d_du;
constraints.C_u   = C_u;
constraints.L     = L;
constraints.u_min = u_min;
constraints.u_max = u_max;

%% ----------------- Salvar dados -----------------
controller_data.A_pos = A_pos;
controller_data.B_pos = B_pos;
controller_data.C_pos = C_pos;
controller_data.nx_pos = nx_pos;

controller_data.A = A; controller_data.B = B; controller_data.C = C;
controller_data.B_v = B_v; controller_data.nx = nx;

controller_data.G = G; controller_data.F = F; controller_data.H = H;
controller_data.Hess = Hess; controller_data.grad_base = grad_base;
controller_data.constraints = constraints;

controller_data.N2 = N2; controller_data.Nu = Nu;
controller_data.ny = ny; controller_data.nu = nu; controller_data.nv = nv;

controller_data.u0 = u0_siso;
controller_data.y0 = y0_siso;

save(save_filename,'controller_data');
fprintf('Arquivo salvo: %s\n\n', save_filename);
end