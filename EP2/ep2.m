#!/usr/bin/octave -qf


function [ind v] = simplex(A, b, c, m, n, x)

    B = zeros(n-m, 1);

    k = 1;
    for j = 1:n
        if (x(j) > 0)
            B(k++) = j;
        endif
    endfor

    [ind v] = simplexR(A, b, c, m, n, x, B);

endfunction


function [ind v] = simplexR(A, b, c, m, n, x, B)
    
    persistent cont = 0;
    ind = 0;
    k = 0;

    B_inv = inv(A(:, B));

    cst_r = zeros(n, 1);

    printf('\nIterando %d\n', cont++);
    printf('-----------\n');

    for j = 1:n
        if (x(j) > 0)
            printf('x%d -> %.5g\n', j, x(j));
        else
            if (size(B_inv) == 0)
                continue;
            endif
            cst_r(j) = c(j) - (transpose(c(B)) * B_inv * A(:, j));
            if (cst_r(j) < 0)
                ind = 1;
                k = j;
            endif
        endif
    endfor

    obj = transpose(c) * x;
    printf('\nValor função objetivo: %d\n', obj);

    printf('\nCustos reduzidos:\n');
    for j = 1:n
        if (x(j) == 0)
            printf('c%d -> %.5g\n', j, cst_r(j));
        endif
    endfor

    if (ind == 0)
        v = x;

    else
        u = B_inv * A(:, k);

        if (u <= 0)
            [ind v] = deal(-1, -u);

        else
            printf('\nSai da base: %d\n', k);

            [theta l] = get_theta(x, u, m, n, B);

            printf('\nTheta*\n');
            printf('%.5g\n', theta);

            printf('\nDireção:\n');
            for j = 1:n-m
                printf('d%d -> %.5g\n', j, -u(j));
            endfor

            printf('\nEntra da base: %d\n', l);

            x(k) = theta;
            for i = 1:n-m
                x(B(i)) -= theta * u(i);
            endfor
            x(l) = 0;

            A(:, B(l)) = A(:, k);
            B(l) = k;
            B = sort(B);

            [ind v] = simplexR(A, b, c, m, n, x, B);

        endif
    endif

endfunction


function [minimo l] = get_theta(x, u, m, n, B)

    [minimo l] = deal(x(B(1)) / u(1), 1);

    for i = 2:n-m
        if (u(i) <= 0) 
            continue;
        endif
        valor = x(B(i)) / u(i);
        if (valor < minimo)
            [minimo l] = deal(valor, i);
        endif
    endfor  

endfunction


function A = le_matriz(arq, m, n)

    A = zeros(m, n);
    for i = 1:m
        for j = 1:n
            A(i, j) = fscanf(arq, '%f', 1);
        endfor
    endfor

endfunction


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

[inv v] = simplex(A, b, c, m, n, x);

if (inv == -1)
    printf('\nO problema é ilimitado!\n');
else
    custo = transpose(c) * v;
    printf('\nSolução ótima encontrada com custo %.5g:\n', custo);
endif

for j = 1:n
    if (inv == -1)
        printf('d');
    else
        printf('x');
    endif
    printf('%d -> %.5g\n', j, v(j));
endfor
