#!/usr/bin/octave --silent

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
    if ! prod(A*x == b)
        erro('Não é uma solução viável. Verifique sua entrada!');
        return
    end
    printf('Simplex: Fase 2\n');
    B = [];
    it = 0; 
    k = 1;
    
    # Monta vetor de índices básicos
    for j = 1:n
        if (x(j) > 0)
            B(k++, 1) = j;
        end
    end
    if (length(B) != m)
        erro('Não é uma solução básica não degenerada. Verifique sua entrada!');
    end
    ind = 1;
    # Toma as colunas de A correspondentes aos índices básicos e calcula B⁻¹
    B_inv = inv(A(:, B));
    
    # Roda o laço enquanto não encontrar solução ou direção ótima
    while ind != 0 && ind != -1
        printf('\nIterando %d\n', it++);
        for i = 1:m
            printf('%d %f\n', B(i), x(B(i)));
        end
        printf('\nValor função objetivo: %f\n\n', transpose(c) * x);
        B_inv = inv(A(:, B));
        # Pré-calcula p a fim de evitar operações desnecessárias
        p = transpose((transpose(c(B))) * B_inv);
        cst_r = zeros(n, 1);
        printf('Custos reduzidos\n');
        for j = 1:n
            if (x(j) == 0)
                # Calcula custo reduzido
                cst_r(j) = c(j) - transpose(p) * A(:, j);
                printf('%d %f\n', j, cst_r(j));
            end
        end
        j = 1;
        while j <= n && cst_r(j) >= 0
            j++;
        end
        if j > n # Se todos cst_r forem positivos, encontramos solução ótima
            printf('\nSolução ótima encontrada com custo %f:\n', transpose(c) * x);
            for j = 1:n
                printf('%d %f\n', j, x(j));
            end
            ind = 0;
            v = x;
        else # cj é negativo, então x não é ótimo; continua o algoritmo
            # Tomamos u como sendo -dB
            u = B_inv * A(:, j);
            if (u <= 0)
                # Se nenhuma componente de u for positiva, então dB > 0.
                # Logo, θ* = +∞ e o custo ótimo será -∞.
                v = zeros(n, 1);
                printf('\nO custo ótimo é -infinito. Direção:\n');
                for i = 1:m
                    printf('%d %f\n', i, u(i));
                    v(B(i)) = u(i);
                end
                ind = -1;
            else
                l = 1;
                theta = inf;
                for i = 1:m
                    if u(i) > 0 && (t = x(B(i)) / u(i)) < theta
                        l = i;
                        theta = t;
                    end
                end
                printf('\nEntra na base: %d\n\nDireção\n', j);
                for i = 1:m
                    if u(i) > 0
                        printf('%d %f\n', B(i), u(i));
                    end
                end
                printf('\nTheta*\n%f\n', theta);
                printf('\nSai da base: %d\n', B(l));
                # Atualiza o valor de x com a nova s.v.b. encontrada
                x(j) = theta;
                x(B(l)) = 0;
                for i = [1:(l - 1), (l + 1):m]
                    x(B(i)) -= theta * u(i);
                end
                # Atualiza vetor de índices básicos, mantendo-o ordenado
                disp(B);
                B(l) = []; 
                B(m, 1) = 0; # apenas para manter o tamanho
                i = m - 1;
                while i > 0 && B(i) > j
                    B(i + 1) = B(i);
                    i--;
                end
                B(i + 1) = j;
                # Atualiza matriz B⁻¹ 
                for i = [1:(l - 1), (l + 1):m]
                    k = -u(i) / u(l);
                    B_inv(i, :) += -k * B_inv(l, :);
                end
                B_inv(l, :) /= u(l);
            end
        end
    end
end

function A = le_matriz(arq, m, n)
    A = zeros(m, n);
    for i = 1:m
        for j = 1:n
            A(i, j) = fscanf(arq, '%f', 1);
        end
    end
end

function erro(msg)
    printf('ERRO: %s\n', msg);
    exit();
end


nome_arq = argv(){1};
arq = fopen(nome_arq, 'r');
m = fscanf(arq, '%f', 1);
n = fscanf(arq, '%f', 1);
A = le_matriz(arq, m, n);
b = le_matriz(arq, m, 1);
c = le_matriz(arq, n, 1);
x = le_matriz(arq, n, 1);

simplex(A, b, c, m, n, x);
