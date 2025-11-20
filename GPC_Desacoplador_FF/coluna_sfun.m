function coluna_sfun(block)
% =========================================================================
% coluna_sfun.m (Level-2 S-Function)
% Modelo dinâmico não linear de uma coluna de destilação binária
% Saídas: [xB; xD]
% =========================================================================

setup(block);

end % coluna_sfun


function setup(block)
    % Ler parâmetro 'p' passado no mask/dialog
    p = block.DialogPrm(1).Data;

    % --- Registrar número de portas ---
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 1;

    % --- Porta de entrada ---
    block.InputPort(1).Dimensions        = 5; % Entradas: L, V, F, zF, qF
    block.InputPort(1).DatatypeID        = 0; % double
    block.InputPort(1).Complexity        = 'Real';
    block.InputPort(1).DirectFeedthrough = false;

    % --- Porta de saída ---
    block.OutputPort(1).Dimensions       = 2; % Somente xB e xD
    block.OutputPort(1).DatatypeID       = 0; % double
    block.OutputPort(1).Complexity       = 'Real';

    % --- Número de parâmetros ---
    block.NumDialogPrms = 1; % Parâmetro 'p'
    
    % --- Estados contínuos ---
    block.NumContStates = 2 * p.NT;

    % --- Tempo de amostragem ---
    block.SampleTimes = [0 0]; % Contínuo

    % --- Métodos ---
    block.RegBlockMethod('InitializeConditions', @InitializeConditions);
    block.RegBlockMethod('Outputs',              @Outputs);
    block.RegBlockMethod('Derivatives',          @Derivatives);
end % setup


function InitializeConditions(block)
    % Inicializar estados
    p  = block.DialogPrm(1).Data;
    NT = p.NT;
    
    % Perfil inicial de composição
    x0_ss = linspace(0.01, 0.99, NT)'; 
    
    % M0 pode ser escalar ou vetor
    if isscalar(p.M0)
        M0_vec = p.M0 * ones(NT,1);
    else
        M0_vec = p.M0(:);
    end

    block.ContStates.Data = [M0_vec .* x0_ss; M0_vec];
end


function Outputs(block)
    % Calcular saídas
    p  = block.DialogPrm(1).Data;
    NT = p.NT;
    
    x_all = block.ContStates.Data;
    Mx = x_all(1:NT);
    M  = x_all(NT+1:end);

    % Composição líquida
    comp_x = Mx ./ (M + eps);
    comp_x = min(max(comp_x, 0), 1);

    xD = comp_x(NT); % Destilado (topo)
    xB = comp_x(1);  % Fundo

    % Somente duas saídas
    block.OutputPort(1).Data = [xD; xB];
end


function Derivatives(block)
    % =========================================================================
    %               FUNÇÃO DERIVATIVES COMPLETA E CORRIGIDA
    % =========================================================================

    % --- Parâmetros e Estados ---
    p = block.DialogPrm(1).Data;
    NT = p.NT;
    NF = p.NF;
    alpha = p.alpha;
    M0 = p.M0;
    tau_h = p.tau_h;

    if isscalar(M0)
        M0_vec = M0 * ones(NT,1);
    else
        M0_vec = M0(:);
    end

    % --- Entradas ---
    u = block.InputPort(1).Data;
    L_in = u(1); V_in = u(2); F = u(3); zF = u(4); qF = u(5);

    % --- Vetor de Estados ---
    x_vec = block.ContStates.Data;
    Mx = x_vec(1:NT);
    M  = x_vec(NT+1:end);

    % --- Composições ---
    comp_x = Mx ./ (M + eps);
    comp_x = min(max(comp_x,0),1);
    comp_y = alpha .* comp_x ./ (1 + (alpha-1).*comp_x);

    % --- [CORREÇÃO 1] Cálculo do Perfil de Vazão de Líquido ---
    L0 = zeros(NT,1);
    L_rect = L_in;                  % Vazão de líquido na seção de enriquecimento (acima da alimentação)
    L_strip = L_in + qF * F;        % Vazão de líquido na seção de esgotamento (abaixo da alimentação)

    % Aplica as vazões às seções corretas da coluna
    % Seção de Esgotamento (dos pratos 2 até o prato de alimentação NF)
    if NF >= 2
        L0(2:NF) = L_strip;
    end
    % Seção de Enriquecimento (dos pratos NF+1 até o topo NT)
    if NF < NT
        L0(NF+1:NT) = L_rect;
    end

    % Adiciona a dinâmica hidráulica
    L = L0 + (M - M0_vec) ./ tau_h;
    L(1) = 0; % Não há líquido entrando no refervedor (prato 1) por cima

    % --- Cálculo do Perfil de Vazão de Vapor (Já estava correto) ---
    V = V_in * ones(NT,1);
    if NF <= NT-1
        V(NF:NT-1) = V(NF:NT-1) + (1-qF)*F;
    end

    % --- [CORREÇÃO 2] Cálculo das Vazões de Produto D e B ---
    % O vapor que chega no condensador (estágio NT) é V(NT-1).
    % O refluxo que sai do condensador é L_in.
    D = V(NT-1) - L_in;
    B = F - D;

    % --- Inicialização das Derivadas ---
    dMx = zeros(NT,1);
    dM  = zeros(NT,1);

    % --- [CORREÇÃO 3] Balanços de Massa ---

    % Fundo (Refervedor - Estágio 1) com retirada de produto B
    dM(1)  = L(2) - V(1) - B;
    dMx(1) = L(2)*comp_x(2) - V(1)*comp_y(1) - B*comp_x(1);

    % Estágios Intermediários
    for i=2:NT-1
        dM(i)  = L(i+1)-L(i) + V(i-1)-V(i);
        dMx(i) = L(i+1)*comp_x(i+1) - L(i)*comp_x(i) + V(i-1)*comp_y(i-1) - V(i)*comp_y(i);
    end

    % Adição da Alimentação no estágio NF
    dM(NF)  = dM(NF)  + F;
    dMx(NF) = dMx(NF) + F*zF;

    % Topo (Condensador - Estágio NT) com retirada de produto D
    % L(NT) é a vazão de refluxo que sai do condensador.
    dM(NT)  = V(NT-1) - L(NT) - D;
    dMx(NT) = V(NT-1)*comp_y(NT-1) - L(NT)*comp_x(NT) - D*comp_x(NT);

    % --- Atribuição Final das Derivadas ---
    block.Derivatives.Data = [dMx; dM];
end