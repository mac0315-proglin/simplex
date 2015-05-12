#!/usr/bin/octave -qf


function [ind v] = simplex(A, b, c, m, n, x)

    ind = cont = 0;

    # Monta vetor de índices básicos
    B = [];
    k = 1;
    for j = 1:n
        if (x(j) > 0)
            B(k++) = j;
        end
    end

    # Valor inicial de B^-1
    B_inv = inv(A(:, B));

    # Custos reduzidos
    cst_r = [];

    # Pré-calcula p para evitar operações desnecessárias
    p = (transpose(c(B)) * B_inv

    while (ind != 0 && ind != -1)

        for j = 1:n
            if (x(j) > 0)
                printf('x%d -> %.5g\n', j, x(j));
            else
                if (size(B_inv) == 0)
                    continue;
                end
                cst_r(j) = c(j) - (p * A(:, j));
                if (cst_r(j) < 0)
                    ind = 1;
                    k = j;
                end
            end
        end

        printf('\nIterando %d\n', cont++);
        printf('-----------\n');

        printf('\nValor função objetivo: %d\n', transpose(c) * x);

        printf('\nCustos reduzidos:\n');
        for j = 1:n
            if (x(j) == 0)
                printf('c%d -> %.5g\n', j, cst_r(j));
            end
        end

        if (ind == 0)
            v = x;

        else
            u = B_inv * A(:, k);

            # Se 
            if (u <= 0)
                [ind v] = deal(-1, -u);

            else
                printf('\nSai da base: %d\n', k);

                [theta l] = calcula_theta(x, u, m, n, B);

                printf('\nTheta*\n');
                printf('%.5g\n', theta);

                printf('\nDireção:\n');
                for j = 1:n-m
                    printf('d%d -> %.5g\n', j, -u(j));
                end

                printf('\nEntra da base: %d\n', l);

                x(k) = theta;
                for i = 1:n-m
                    x(B(i)) -= theta * u(i);
                end
                x(l) = 0;

                A(:, B(l)) = A(:, k);

                B_inv = calcula_inv(B_inv, u, l)
            end
        end
    end
end


function [theta l] = calcula_theta(x, u, m, n, B)

    theta = inf;

    for i = 1:m
        if (u(i) <= 0) 
            continue;
        end
        valor = x(B(i)) / u(i);
        if (valor < theta)
            [theta l] = deal(valor, i);
        end
    end
end


function B_inv = calcula_inv(B_inv, u, l)

    for j = [1:(l-1), (l+1):m]
        k = -u(i)/u(l)
        B_inv(i) += -k * B_inv(l)
    end
    B_inv(l) /= u(l)
end


function A = le_matriz(arq, m, n)

    A = zeros(m, n);
    for i = 1:m
        for j = 1:n
            A(i, j) = fscanf(arq, '%f', 1);
        end
    end
end


# ------------------------------------------------------------- #

printf('===========================\n')
printf('===   Simplex: Fase 2   ===\n')
printf('===========================\n');

nome_arq = argv(){1};
arq = fopen(nome_arq, 'r');

m = fscanf(arq, '%f', 1)
n = fscanf(arq, '%f', 1)

A = le_matriz(arq, m, n)
b = le_matriz(arq, m, 1)
c = le_matriz(arq, n, 1)
x = le_matriz(arq, n, 1)

if (not prod(A*x == b) or )
    printf('x NÃO É UMA SOLUÇÃO VIÁVEL! Verifique sua entrada.');
    return
end

[inv v] = simplex(A, b, c, m, n, x);

if (inv == -1)
    printf('\nO problema tem custo ótimo -∞\n');
else
    custo = transpose(c) * v;
    printf('\nSolução ótima encontrada com custo %.5g:\n', custo);
end

for j = 1:n
    if (inv == -1)
        printf('d');
    else
        printf('x');
    end
    printf('%d -> %.5g\n', j, v(j));
end
