function [L, D] = rozklad_cholesky_LDLh_optimized(A)
% rozklad_cholesky_LDLh_optimized Optymalny rozkład LDL^h dla macierzy pięciodiagonalnej hermitowskiej i dodatnio określonej.
%
% Dane wejściowe:
%   - A: macierz pięciodiagonalna hermitowska i dodatnio określona.
%
% Dane wyjściowe:
%   - L: macierz dolnotrójkątna, gdzie główna przekątna to jedynki.
%   - D: macierz diagonalna, zawierająca wartości D(k).
%
% Funkcja wykorzystuje algorytm rozkładu LDL^h zoptymalizowany dla macierzy pięciodiagonalnych.

    if ~ifProperMatrix(A) 
        disp("Macierz nie spełnia założeń")
        return
    end
    A_diag = diag(A, 0);
    A_lower1 = diag(A, -1);
    A_lower2 = diag(A, -2);
    n = length(A_diag); 
    L_lower1 = zeros(n-1, 1);  
    L_lower2 = zeros(n-2, 1);  
    D = zeros(n, 1);  
    % Rozkład LDL^h
    for k = 1:n
        if k == 1
            D(k) = A_diag(k);
        elseif k == 2
            D(k) = A_diag(k) - (abs(L_lower1(k-1))^2 * D(k-1));
        else
            D(k) = A_diag(k) - (abs(L_lower1(k-1)))^2 * D(k-1) - (abs(L_lower2(k-2))^2 * D(k-2));
        end
        if k < n 
            if k > 1
                L_lower1(k) = (A_lower1(k) - (conj(L_lower1(k-1))* D(k-1)*L_lower2(k-1))) / D(k);  
            else
                 L_lower1(k) = A_lower1(k) / D(k);
            end
        end
        if k < n-1
            L_lower2(k) = A_lower2(k) / D(k);            
        end
    end

    L = eye(n); 
    L(sub2ind([n, n], 2:n, 1:n-1)) = L_lower1; 
    L(sub2ind([n, n], 3:n, 1:n-2)) = L_lower2;  
    D = diag(D);
end