#!/usr/bin/octave --silent

%
%  Método simplex revisado.
%
%  Recebe como parâmetros uma matriz A ∈ R^(m x n), vetores b ∈ R^m e c ∈ R^n
%  representando o seguinte problema no formato padrão:
%
%      Minimizar c'x
%      Sujeito a Ax = b
%                 x ≥ 0
%
%  A função aplica o método simplex (implementação revisada) e devolve como
%  resultado um indicador 'ind', bem como possíveis valores 'x' e 'd':
%    - Caso o problema tenha solução finita, ind == 0 e x será a solução viável
%      básica que minimiza o custo c'x;
%    - Caso seja ilimitado, ind == -1 e d será a direção viável na qual
%      o problema tende ao custo -∞;
%    - Por fim, se o problema for inviável, ind == 1.
%
function [ind x d] = simplex(A, b, c, m, n)

    % Multiplicando por -1 as restrições em que b(i) < 0
    for i = 1:m
        if b(i) < 0
            A(i, :) *= -1;
            b(i) *= -1;
        end
    end

    % Constrói o problema auxiliar
    A_aux = [A, eye(m)];
    c_aux = [zeros(n, 1); ones(m, 1)];
    x = [zeros(n, 1); b];

    printf('\n===========================')
    printf('\n===   Simplex: Fase 1   ===')
    printf('\n===========================');;
    [ind x u B B_inv] = simplex_body(A_aux, b, c_aux, m, n + m, x);

    % O problema da fase 1 sempre é viável (pois há uma solução trivial)
    % e sempre tem custo ótimo finito (pois é a soma de variáveis
    % não-negativas); se o custo ótimo da fase 1 for estritamente
    % positivo, então o problema original é inviável
    if transpose(c_aux) * x > 0
        ind = 1;

    else
        % Retira da base as variáveis artificiais
        for l = 1:m

            % (l-ésima linha de B⁻¹) A
            v = B_inv(l, :) * A;

            for j = 1:m
                if v(j) == 0
                    break;
                end
            end

            if j > m
                % v é nulo, então a l-ésima restrição é redundante e é removida
                A(l, :) = [];
                B(l) = [];
                m -= 1;

                printf('A restrição #%d do problema original era redundante.\n', l);
                
                % TODO: atualizar B⁻¹ eficientemente
                B_inv = inverse(A(:, B));

            else 
                % v(j) não é nulo, então x_j entra na base, x_l sai
                u = B_inv * A(:, j);
                for i = [1:(l-1), (l+1):m]
                    r = -u(i) / u(l);
                    B_inv(i, :) += r * B_inv(l, :);
                end
                B_inv(l, :) /= u(l);
                % Atualiza vetor de índices básicos
                B(l) = j;
            end
        end

        x = x(1:n);

        printf('\n===========================')
        printf('\n===   Simplex: Fase 2   ===')
        printf('\n===========================');
        [ind x u B] = simplex_body(A, b, c, m, n, x);

        if ind == -1
            % Monta vetor de direção
            d = zeros(n, 1);
            for i = 1:m
                d(B(i)) = -u(i);
            end
        end
    end
end

%
% Núcleo das iterações das fases do método simplex.
%
function [ind x u B B_inv] = simplex_body(A, b, c, m, n, x)

    ind = -2;
    cont = k = 0;
    B = cst_r = d = u = [];

    % Monta vetor B de índices básicos
    for j = 1:n
        if x(j) > 0
            B(++k) = j;
        end
    end

    % Toma as colunas de A correspondentes aos índices básicos e calcula B⁻¹
    B_inv = inv(A(:, B));

    % Roda o laço enquanto não encontrar solução ou direção ótima
    while (ind ~= 0) && (ind ~= -1)

        printf('\n\n-------------');
        printf('\n- Iterando %d', cont++);
        printf('\n-------------\n');
        printf('\n> Valor da função objetivo: %.5g\n', transpose(c) * x);

        % Pré-calcula p a fim de evitar operações desnecessárias
        p = (transpose(c(B))) * B_inv;
      
        k = 0;
        for j = 1:n
            if any(B(:) == j) % se está na base
                printf('x%d -> %.5g\n', j, x(j));
            else
                % Calcula custo reduzido
                cst_r(j) = c(j) - (p * A(:, j));

                % Se negativo, x não é ótimo; continua o algoritmo
                % usando a regra do menor índice
                if cst_r(j) < 0 && k == 0
                    k = j;
                end
            end
        end

        printf('\n> Custos reduzidos:\n');
        for j = 1:n
            if ~ any(B(:) == j)  % se não está na base
                printf('c%d -> %.5g\n', j, cst_r(j));
            end
        end

        if k == 0
            % Se todos cst_r forem não negativos, encontramos solução ótima
            ind = 0;

        else
            % Tomamos u como sendo -dB
            u = B_inv * A(:, k);

            if u <= 0
                % Se nenhuma componente de u for positiva, então dB > 0.
                % Logo, θ* = +∞ e o custo ótimo será -∞.
                ind = -1;

            else
                [theta l] = calcula_theta(x, u, m, n, B);
                printf('\n> Theta*: (%.5g)\n', theta);

                printf('\n> Direção:\n');
                for i = 1:m
                    printf('d%d -> %.5g\n', B(i), -u(i));
                end

                printf('\n> Sai da base: (%d)', B(l));
                printf('\n> Entra na base: (%d)\n', k);

                % Atualiza o valor de x com a nova s.v.b encontrada
                x(k) = theta;
                for i = 1:m
                    x(B(i)) -= theta * u(i);
                end

                % Calcula eficientemente a nova B⁻¹ utilizando de operações
                % elementares entre as linhas e o vetor u (simplex revisado).
                for i = [1:(l-1), (l+1):m]
                    r = -u(i) / u(l);
                    B_inv(i, :) += r * B_inv(l, :);
                end
                B_inv(l, :) /= u(l);

                % Atualiza vetor de índices básicos
                B(l) = k;
            end
        end
    end
end

%
%  Recebe um ponto x em R^n, direção u e um conjunto B de m índices.
%  Devolve θ*, isto é, o menor valor da razão x(B(i))/u(i), i = 1,...,m.
%
function [theta l] = calcula_theta(x, u, m, n, B)

    theta = inf;

    for i = 1:m
        if u(i) > 0
            valor = x(B(i)) / u(i);
            if valor < theta
                [theta l] = deal(valor, i);
            end
        end
    end
end

%
%  Lê arquivo e devolve matriz de tamanho m x n.
%
function A = le_matriz(arq, m, n)

    A = zeros(m, n);
    for i = 1:m
        for j = 1:n
            A(i, j) = fscanf(arq, '%f', 1);
        end
    end
end

%
%  Imprime mensagem de erro.
%
function erro(msg)

    printf('\n[ERRO]\n> %s', msg);
    exit();
end

% -----------------------------------------------------------------------------

% Abre o arquivo
nome_arq = argv(){1};
arq = fopen(nome_arq, 'r');

% Leitura da entrada
m = fscanf(arq, '%f', 1);
n = fscanf(arq, '%f', 1);
A = le_matriz(arq, m, n);
b = le_matriz(arq, m, 1);
c = le_matriz(arq, n, 1);
fclose(arq);

% Chamada da função
[ind v] = simplex(A, b, c, m, n);

printf('\n\n-------------');
printf('\n- RESULTADO -');
printf('\n-------------\n');

switch (ind)
    case 0
        printf('\n> Solução ótima encontrada com custo %.5g:\n',
                    transpose(c) * v);
    case -1
        printf('\n> O problema tem custo ótimo -∞\n');
        printf('\n> Direção viável geradora:\n');
    case 1
        printf('\n> O problema é inviável!\n');
end

% Imprime resultado, isto é, a solução (ou direção) ótima
if ind ~= 1
    for j = 1:n
        if ind == -1
            printf('d');
        else
            printf('x');
        end
        printf('%d -> %.5g\n', j, v(j));
    end
end