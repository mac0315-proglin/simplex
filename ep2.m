#!/usr/bin/octave -qf

#
#  Método simplex revisado.
#
#  Recebe como parâmetros uma matriz A ∈ R^(m x n), vetores b ∈ R^m e c ∈ R^n
#  e um ponto x ∈ R^n, representando o seguinte problema no formato padrão:
#
#      Minimizar c'x
#      Sujeito a Ax = b
#                 x ≥ 0
#
#  A função aplica o método simplex (implementação revisada) e devolve como
#  resultado um indicador 'ind' e um valor 'v':
#    - Caso o problema tenha solução finita, ind == 0 e v será a solução viável
#      básica que minimiza o custo c'x.
#    - Caso contrário, ind == -1 e v será a direção viável d na qual o problema
#      tende ao custo -∞.
#
function [ind v] = simplex(A, b, c, m, n, x)

    cont = k = 0;
    B = cst_r = [];

    # Monta vetor de índices básicos
    for j = 1:n
        if x(j) > 0
            B(++k) = j;
        end
    end

    if k > m
        erro('x não é solução viável básica\n');
    end

    # Toma as colunas de A correspondentes aos índices básicos e calcula B⁻¹
    B_inv = inv(A(:, B));

    # Roda o laço enquanto não encontrar solução ou direção ótima
    while (ind ~= 0) and (ind ~= -1)

        printf('\nIterando %d\n', cont++);
        printf('-----------\n');
        printf('\nValor função objetivo: %f\n', transpose(c) * x);

        # Pré-calcula p a fim de evitar operações desnecessárias
        p = (transpose(c(B))) * B_inv;

        ind = 0;
        for j = 1:n
            if x(j) > 0
                printf('x%d -> %.5g\n', j, x(j));
            else
                # VERIFICAR se esta condição é necessária (caso extremo)
                #if (size(B_inv) == 0)
                #    continue;
                #end
                # Calcula custo reduzido
                cst_r(j) = c(j) - (p * A(:, j));

                # Se negativo, x não é ótimo; continua o algoritmo
                if cst_r(j) < 0
                    ind = 1;
                    k = j;
                end
            end
        end

        printf('\nCustos reduzidos:\n');
        for i = 1:m
            if x(B(i)) == 0
                printf('c%d -> %.5g\n', i, cst_r(i));
            end
        end

        if ind == 0
            # Se todos cst_r forem positivos, encontramos solução ótima
            v = x;

        else
            # Tomamos u como sendo -dB
            u = B_inv * A(:, k);

            if u <= 0
                # Se nenhuma componente de u for positiva, então dB > 0.
                # Logo, θ* = +∞ e o custo ótimo será -∞.
                [ind v] = deal(-1, -u);

            else
                printf('\nSai da base: %d\n', k);

                [theta l] = calcula_theta(x, u, m, n, B);

                if theta == 0
                    erro('Encontrada s.v.b degenerada. Verifique entrada.\n');
                end

                printf('\nTheta*\n');
                printf('%.5g\n', theta);

                printf('\nDireção:\n');
                for i = 1:m
                    printf('d%d -> %.5g\n', i, -u(i));
                end

                printf('\nEntra da base: %d\n', l);

                # Atualiza o valor de x com a nova s.v.b encontrada
                x(k) = theta;
                for i = 1:m
                    x(B(i)) -= theta * u(i);
                end

                # Atualiza matrizes A e B⁻¹
                A(:, B(l)) = A(:, k);
                B_inv = calcula_inv(B_inv, u, m, l);

                # Atualiza e ordena vetor de índices básicos
                B(:, B(l)) = A(:, k);
                sort(B);
            end
        end
    end
end

#
#  Recebe um ponto x em R^n, direção u e um conjunto B de m índices.
#  Devolve θ*, isto é, o menor valor da razão x(B(i))/u(i), i = 1,...,m.
#
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

#
#  Recebe uma matriz B⁻¹, um vetor u e índices m e l.
#  Calcula a nova inversa de B (com a coluna de índice l sendo agora u) por
#  meio de operações elementares entre as linhas (simplex revisado).
#  Devolve matriz B⁻¹ atualizada.
#
function B_inv = calcula_inv(B_inv, u, m, l)

    for i = [1:(l-1), (l+1):m]
        k = -u(i) / u(l);
        B_inv(i) += -k * B_inv(l);
    end
    
    B_inv(l) /= u(l);
end

#
#  Lê arquivo e devolve matriz de tamanho m x n.
#
function A = le_matriz(arq, m, n)

    A = zeros(m, n);
    for i = 1:m
        for j = 1:n
            A(i, j) = fscanf(arq, '%f', 1);
        end
    end
end

#
#  Imprime mensagem de erro.
#
function erro(msg)

    puts('ERRO: ', msg);
    exit();
end

# -----------------------------------------------------------------------------

printf('===========================\n')
printf('===   Simplex: Fase 2   ===\n')
printf('===========================\n');

# Abre o arquivo
nome_arq = argv(){1};
arq = fopen(nome_arq, 'r');

# Leitura dos valores
m = fscanf(arq, '%f', 1)
n = fscanf(arq, '%f', 1)
A = le_matriz(arq, m, n)
b = le_matriz(arq, m, 1)
c = le_matriz(arq, n, 1)
x = le_matriz(arq, n, 1)

if not (A*x == b)
    erro('x não é solução viável. Verifique sua entrada.');
end

# Chamada da função
[ind v] = simplex(A, b, c, m, n, x);

if ind == -1
    printf('\nO problema tem custo ótimo -∞');
    printf('\nDireção viável geradora:\n');
else
    printf('\nSolução ótima encontrada com custo %.5g:\n', transpose(c) * v);
end

# Imprime resultado, isto é, a solução (ou direção) ótima
for j = 1:n
    if ind == -1
        printf('d');
    else
        printf('x');
    end
    printf('%d -> %.5g\n', j, v(j));
end
