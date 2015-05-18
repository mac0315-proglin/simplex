#!/usr/bin/octave --silent

%
%  Método simplex revisado.
%
%  Recebe como parâmetros uma matriz A ∈ R^(m x n), vetores b ∈ R^m e c ∈ R^n
%  e um ponto x ∈ R^n, representando o seguinte problema no formato padrão:
%
%      Minimizar c'x
%      Sujeito a Ax = b
%                 x ≥ 0
%
%  A função aplica o método simplex (implementação revisada) e devolve como
%  resultado um indicador 'ind' e um valor 'v':
%    - Caso o problema tenha solução finita, ind == 0 e v será a solução viável
%      básica que minimiza o custo c'x.
%    - Caso contrário, ind == -1 e v será a direção viável d na qual o problema
%      tende ao custo -∞.
%
function [ind v] = simplex(A, b, c, m, n, x)

    ind = 1;
    cont = k = 0;
    B = cst_r = [];

    % Monta vetor B de índices básicos
    for j = 1:n
        if x(j) > 0
            B(++k) = j;
        end
    end

    if k > m
        erro('x não é solução viável básica\n');
    end

    % Toma as colunas de A correspondentes aos índices básicos e calcula B⁻¹
    B_inv = inv(A(:, B));

    % Roda o laço enquanto não encontrar solução ou direção ótima
    while (ind ~= 0) && (ind ~= -1)

        printf('\n\n-------------');
        printf('\n- Iterando %d', cont++);
        printf('\n-------------\n');
        printf('\n> Valor função objetivo: %.5g\n', transpose(c) * x);

        % Pré-calcula p a fim de evitar operações desnecessárias
        p = (transpose(c(B))) * B_inv;

        ind = num_positivos = 0;
        for j = 1:n
            if x(j) > 0
                printf('x%d -> %.5g\n', j, x(j));
                num_positivos++;
            else
                % Calcula custo reduzido
                cst_r(j) = c(j) - (p * A(:, j));

                % Se negativo, x não é ótimo; continua o algoritmo
                if cst_r(j) < 0
                    ind = 1;
                    k = j;
                end
                
            end
        end
        if num_positivos < m
            erro('Encontrada s.v.b degenerada. Verifique sua entrada.');
        end

        printf('\n> Custos reduzidos:\n');
        for j = 1:n
            if x(j) == 0
                printf('c%d -> %.5g\n', j, cst_r(j));
            end
        end

        if ind == 0
            % Se todos cst_r forem positivos, encontramos solução ótima
            v = x;

        else
            % Tomamos u como sendo -dB
            u = B_inv * A(:, k);

            if u <= 0
                % Se nenhuma componente de u for positiva, então dB > 0.
                % Logo, θ* = +∞ e o custo ótimo será -∞.

                % Monta vetor de direção
                d = zeros(n, 1);
                for i = 1:m
                    d(B(i)) = -u(i);
                end
                d(k) = 1;

                [ind v] = deal(-1, d);

            else
                [theta l] = calcula_theta(x, u, m, n, B);
                printf('\n> Theta*: (%.5g)\n', theta);

                printf('\n> Direção:\n');
                for i = 1:m
                    printf('d%d -> %.5g\n', i, -u(i));
                end

                printf('\n> Sai da base: (%d)', k);
                printf('\n> Entra da base: (%d)\n', l);

                % Atualiza o valor de x com a nova s.v.b encontrada
                x(k) = theta;
                for i = 1:m
                    x(B(i)) -= theta * u(i);
                end

                % Calcula eficientemente a nova B⁻¹ utilizando de operações
                % elementares entre as linhas e o vetor u (simplex revisado).
                for i = [1:(l-1), (l+1):m]
                    r = -u(i) / u(l);
                    B_inv(i) += -r * B_inv(l);
                end
                B_inv(l) /= u(l);

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

printf('\n===========================')
printf('\n===   Simplex: Fase 2   ===')
printf('\n===========================');

% Abre o arquivo
nome_arq = argv(){1};
arq = fopen(nome_arq, 'r');

% Leitura da entrada
m = fscanf(arq, '%f', 1);
n = fscanf(arq, '%f', 1);
A = le_matriz(arq, m, n);
b = le_matriz(arq, m, 1);
c = le_matriz(arq, n, 1);
x = le_matriz(arq, n, 1);

if not (A*x == b)
    erro('x não é solução viável. Verifique sua entrada.');
end

% Chamada da função
[ind v] = simplex(A, b, c, m, n, x);

printf('\n\n-------------');
printf('\n- RESULTADO -');
printf('\n-------------\n');

if ind == -1
    printf('\n> O problema tem custo ótimo -∞\n');
    printf('\n> Direção viável geradora:\n');
else
    printf('\n> Solução ótima encontrada com custo %.5g:\n', transpose(c) * v);
end

% Imprime resultado, isto é, a solução (ou direção) ótima
for j = 1:n
    if ind == -1
        printf('d');
    else
        printf('x');
    end
    printf('%d -> %.5g\n', j, v(j));
end
