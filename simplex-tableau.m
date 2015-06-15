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
    B = (n + 1):(n + m);
    B_inv = eye(m);

    printf('\n===========================')
    printf('\n===   Simplex: Fase 1   ===')
    printf('\n===========================');;
    [ind x u B B_inv] = simplex_body(A_aux, b, c_aux, m, n + m, x, B, B_inv);

    % O problema da fase 1 sempre é viável (pois há uma solução trivial)
    % e sempre tem custo ótimo finito (pois é a soma de variáveis
    % não-negativas); se o custo ótimo da fase 1 for estritamente
    % positivo, então o problema original é inviável
    if transpose(c_aux) * x > 0
        ind = 1;
    else
        % Retira da base as variáveis artificiais
        l = 1;
        while l <= m
            if B(l) <= n % não é artificial, nada a fazer
                l++;
                continue;
            end
            imprime_tableau(A_aux, B, b, B_inv, c_aux);
            % v = (l-ésima linha de B⁻¹)  A
            v = B_inv(l, :) * A;

            % procurar elemento não nulo de v
            j = 1;
            while j <= n && v(j) == 0
                j++;
            end

            if j > n
                % v é nulo, então a l-ésima restrição é redundante e é removida
                A(l, :) = [];
                A_aux(l, :) = [];
                b(l) = [];
                printf('\nRemovida restrição redundante do problema original.\n', l);
                B(l) = [];
                m -= 1;
                B_inv(l, :) = [];
                B_inv(:, l) = [];
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
                l++;
            end
        end
        
        imprime_tableau(A_aux, B, b, B_inv, c_aux);

        x = x(1:n);

        printf('\n===========================')
        printf('\n===   Simplex: Fase 2   ===')
        printf('\n===========================');
        [ind x d B] = simplex_body(A, b, c, m, n, x, B, B_inv);

    end
end

%
% Núcleo das iterações das fases do método simplex.
%
function [ind x d B B_inv] = simplex_body(A, b, c, m, n, x, B, B_inv)
    ind = -2;
    cont = k = 0;
    cst_r = d = u = [];

    % Roda o laço enquanto não encontrar solução ou direção ótima
    while (ind ~= 0) && (ind ~= -1)

        printf('\n\n-------------');
        printf('\n- Iterando %d', cont++);
        printf('\n-------------\n');
        printf('\n> Valor da função objetivo: %.5g\n', transpose(c) * x);

        imprime_tableau(A, B, b, B_inv, c);
        % Pré-calcula p a fim de evitar operações desnecessárias
        p = (transpose(c(B))) * B_inv;
        k = 0;
        for j = 1:n
            if any(B(:) == j) % se está na base
                printf('x%s -> %.5g\n', subscrito(j), x(j));
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
                printf('C%s -> %.5g\n', subscrito(j), cst_r(j));
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
                    printf('d%s -> %.5g\n', subscrito(B(i)), -u(i));
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

function [a] = algs(i) 
    switch(i)
        case 1
            a = "₁";
        case 2
            a = "₂";
        case 3
            a = "₃";
        case 4
            a = "₄";
        case 5
            a = "₅";
        case 6
            a = "₆";
        case 7
            a = "₇";
        case 8
            a = "₈";
        case 9
            a = "₉";
        otherwise
            a = "₀";
    end
end

function [subs] = subscrito_r(i)
    algarismo = mod(i, 10);
    if i < 10
        subs = algs(algarismo);
    else
        subs = [subscrito_r(idivide(i, 10)), algs(algarismo)];
    end
end

function [subs] = subscrito(i)
    subs = subscrito_r(i);
    if i < 10
        subs = [subs, " "];
    end
end

function imprime_tableau(A, B, b, B_inv, c)
    [m n] = size(A);
    hifens = ["      ", repmat("-", [1, (8 * n + 16)]), "\n"];
    B_inv_b = B_inv * b;
    B_inv_A = B_inv * A;
    printf("\n\n                    ");
    for j=1:n
        printf("   x%s  ", subscrito(j));
    end
    printf("\n%s", hifens);
    printf("      | %7.3g  |  ", -transpose(c(B)) * B_inv_b);
    printf("%7.3g ", transpose(c) - transpose(c(B)) * B_inv_A);
    printf(" |\n");
    printf("%s", hifens);
    for i=1:m
        printf("x%s = | %7.3g  |  ", subscrito(B(i)), B_inv_b(i));
        printf("%7.3g ", B_inv_A(i, :));
        printf(" |\n");
    end
    printf("%s\n", hifens);
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
            printf('d%s -> %.5g\n', subscrito(j), d(j));
        else
            printf('x%s -> %.5g\n', subscrito(j), x(j));
        end
    end
end

