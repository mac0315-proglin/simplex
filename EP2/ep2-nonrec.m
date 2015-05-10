#!/usr/bin/octave --silent

function [ind v] = simplex(A, b, c, m, n, x)
    printf('Simplex: Fase 2\n\n');
    # // B = zeros(n - m, 1);
    B = [];
    it = 0; 
    k = 1;
    for j = 1:n
        if (x(j) > 0)
            B(k++, 1) = j;
        end
    end
    if (length(B) < n - m)
        printf("NÃO É UMA SOLUÇÃO VIÁVEL BÁSICA! VERIFIQUE SUA ENTRADA!");
        return;
    end
    while true
        printf("Iterando %d\n", it++);
        for i = 1:m
            printf("%d %f\n", B(i), x(B(i)));
        end
        printf('\nValor função objetivo: %f\n\n', transpose(c) * x);
        B_inv = inv(A(:, B));
        cst_r = zeros(n, 1);
        printf('Custos reduzidos\n');
        for j = 1:n
            if (x(j) == 0)
                cst_r(j) = c(j) - transpose(c(B)) * B_inv * A(:, j);
                printf('%d %f\n', j, cst_r(j));
            end
        end
        if (min(cst_r) >= 0)
            printf('\nSolução ótima encontrada com custo %f:\n', transpose(c) * x);
            for j = 1:n
                printf('%d %f\n', j, x(j));
            end
            return
        end
        j = 1
        while cst_r(j) >= 0
            j++;
        end
        u = B_inv * A(:, j);
        if (max(u) <= 0)
            printf('\nO custo ótimo é -infinito. Direção:\n');
            for j = 1:n
                printf('%d %f\n', j, u(j));
            end
            disp(B_inv);
            return
        else
            l = 1;
            theta = inf
            for i = 1:m
                if u(i) > 0 && (t = x(B(i)) / u(i)) < theta
                    l = i;
                    theta = t;
                end
            end
            printf("\nEntra na base: %d\n\nDireção\n", j);
            for i = 1:m
                if u(i) > 0
                    printf("%d %f\n", B(i), u(i));
                end
            end
            printf("\nTheta*\n%f\n", theta);
            printf("\nSai da base: %d\n", B(l));
            x(j) = theta;
            x(B(l)) = 0;
            for i = [1:(l - 1), (l + 1):m]
                x(B(i)) -= theta * u(i);
            end
            B(l) = []; # // removendo elemento do vetor
            B(m, 1) = j;
            input("Digite ENTER:");
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


nome_arq = argv(){1};
arq = fopen(nome_arq, 'r');
m = fscanf(arq, '%f', 1);
n = fscanf(arq, '%f', 1);
A = le_matriz(arq, m, n);
b = le_matriz(arq, m, 1);
c = le_matriz(arq, n, 1);
x = le_matriz(arq, n, 1);

simplex(A, b, c, m, n, x);
