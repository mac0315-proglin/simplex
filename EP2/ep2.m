#!/usr/bin/octave -qf


function [ind v] = simplex(A, b, c, m, n, x)

    cont = 0;
    B = cst_r = [];

    # Monta vetor de índices básicos
    k = 0;
    for j = 1:n
        if x(j) > 0
            B(++k) = j;
        end
    end

    if k ~= m
        printf('\nERRO: x não é solução básica\n');
        exit();
    end

    # Valor inicial de B^-1
    B_inv = inv(A(:, B));

    # TODO: existe um jeito de deixar a condição mais bonita?
    while (ind ~= 0) and (ind ~= -1)

        printf('\nIterando %d\n', cont++);
        printf('-----------\n');
        printf('\nValor função objetivo: %f\n', transpose(c) * x);

        # Pré-calculo a fim de evitar operações desnecessárias
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

                # Se negativo, continua o algoritmo
                if cst_r(j) < 0
                    ind = 1;
                    k = j;
                end
            end
        end

        printf('\nCustos reduzidos:\n');
        for i = 1:m
            if (x(B(i)) == 0)
                printf('c%d -> %.5g\n', i, cst_r(i));
            end
        end

        if ind == 0
            # Se todos cst_r forem positivos, encontramos uma solução ótima
            v = x;

        else
            # Tomamos u como sendo -dB
            u = B_inv * A(:, k);

            # Se nenhuma componente de u for positiva, então dB > 0.
            # Logo, θ* = +∞ e o custo ótimo será -∞.
            if u <= 0
                [ind v] = deal(-1, -u);

            else
                printf('\nSai da base: %d\n', k);

                [theta l] = calcula_theta(x, u, m, n, B);

                if theta == 0
                    printf('\nERRO: Encontrada solução básica degenerada. ');
                    printf('Verifique sua entrada.\n');
                    exit();
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

                # Atualiza matrizes A e B^-1
                A(:, B(l)) = A(:, k);
                B_inv = calcula_inv(B_inv, u, m, l);
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
        if valor < theta
            [theta l] = deal(valor, i);
        end
    end
end


function B_inv = calcula_inv(B_inv, u, m, l)

    for i = [1:(l-1), (l+1):m]
        k = -u(i) / u(l);
        B_inv(i) += -k * B_inv(l);
    end
    
    B_inv(l) /= u(l);

end


function A = le_matriz(arq, m, n)

    A = zeros(m, n);
    for i = 1:m
        for j = 1:n
            A(i, j) = fscanf(arq, '%f', 1);
        end
    end
end


# --------------------------------------------------------------------------- #

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

if (not (A*x == b))
    printf('ERRO: x não é solução viável. Verifique sua entrada.');
    exit()
end

[inv v] = simplex(A, b, c, m, n, x);

if inv == -1
    printf('\nO problema tem custo ótimo -∞\n');
else
    printf('\nSolução ótima encontrada com custo %.5g:\n', transpose(c) * v);
end

for j = 1:n
    if inv == -1
        printf('d');
    else
        printf('x');
    end
    printf('%d -> %.5g\n', j, v(j));
end
