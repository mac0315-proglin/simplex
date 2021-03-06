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
    
    d = [];

    % Multiplicamos por -1 as restrições em que b(i) < 0.
    % Assim, garante-se que b >= 0.
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
    B = (n + 1):(n + m);
    B_inv = eye(m);

    % Aplica a Fase 1 no problema auxiliar
    printf('\n===========================')
    printf('\n===   Simplex: Fase 1   ===')
    printf('\n===========================');;
    [ind x u B B_inv] = simplex_body(A_aux, b, c_aux, m, n + m, x, B, B_inv, n);

    if transpose(c_aux) * x > 0
        % Se o custo ótimo da Fase 1 for positivo,
        % então o problema original é inviável.
        ind = 1;

    else
        % Caso contrário, retiramos da base as variáveis artificiais

        % Percorre-se os índices da base
        l = 1;
        while l <= m

            % Se não for artificial, nada a fazer
            if B(l) <= n
                l++;
                continue;
            end

            % v = l-ésima coluna de B⁻¹A
            v = B_inv(l, :) * A;

            % Procura-se elemento não nulo de v
            j = 1;
            while j <= n && v(j) == 0
                j++;
            end

            if j > n
                % v é nulo, então a l-ésima restrição
                % é redundante e pode ser removida.
                A(l, :) = [];
                B(l) = [];
                B_inv(l, :) = [];
                B_inv(:, l) = [];
                m -= 1;

                printf('\nRemovida restrição redundante.\n', l);                

            else 
                % v(j) não é nulo, então x(j) entra na base, x(l) sai

                u = B_inv * A(:, j);

                % Atualiza B⁻¹ pelo método revisado
                for i = [1:(l-1), (l+1):m]
                    r = -u(i) / u(l);
                    B_inv(i, :) += r * B_inv(l, :);
                end
                B_inv(l, :) /= u(l);

                B(l++) = j;
            end
        end

        % Elimina as variáveis artificiais do jogo
        x = x(1:n);

        % Aplica a Fase 2 do método simplex na base e s.v.b encontradas
        printf('\n===========================')
        printf('\n===   Simplex: Fase 2   ===')
        printf('\n===========================');
        [ind x d B] = simplex_body(A, b, c, m, n, x, B, B_inv, n);
    end
end

%
% Núcleo das iterações das fases do método simplex.
%
function [ind x d B B_inv] = simplex_body(A, b, c, m, n, x, B, B_inv, max_pivot)

    ind = -2;
    cont = k = 0;
    cst_r = d = u = [];

    % Roda o laço enquanto não encontrar solução ou direção ótima
    while (ind ~= 0) && (ind ~= -1)

        printf('\n-------------');
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

        if k == 0 || k > max_pivot
            % Se todos cst_r forem não negativos, encontramos solução ótima
            ind = 0;

        else
            % Tomamos u como sendo -dB
            u = B_inv * A(:, k);

            if u <= 0
                % Se nenhuma componente de u for positiva, então dB > 0.
                % Logo, θ* = +∞ e o custo ótimo será -∞.
                ind = -1;

                % Monta vetor de direção
                d = zeros(n, 1);
                for i = 1:m
                    d(B(i)) = -u(i);
                end
                d(k) = 1;

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
[ind x d] = simplex(A, b, c, m, n);

printf('\n\n-------------');
printf('\n- RESULTADO -');
printf('\n-------------\n');

switch (ind)
    case 0
        printf('\n> Solução ótima encontrada com custo %.5g:\n',
                    transpose(c) * x);
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
            printf('d%d -> %.5g\n', j, d(j));
        else
            printf('x%d -> %.5g\n', j, x(j));
        end
    end
end
