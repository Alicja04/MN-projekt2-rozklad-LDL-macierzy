function [L, D] = matrix_decomposition(A)
%matrix_decomposition Wyznaczenie rozkładu LDL^h dla macierzy hermitowskiej, dodatnio określonej.
%                     Ponadto zakłada się pięciodiagonalność, jednak algorytm można zastosować również do macierzy
%                     niespełniających tego warunku.
%
% Dane wejściowe:
%   - A: macierz hermitowska i dodatnio określona (wymuszany jest również warunek pięciodiagonalności).
%
% Dane wyjściowe:
%   - L: macierz dolnotrójkątna.
%   - D: macierz diagonalna, zawierająca wartości D(k).
%
% Funkcja wyznacza rozkład LDL^h przy użyciu standardowego algorytmu.
    
    [n, ~] = size(A);
    
    if ~ifProperMatrix(A) 
        disp("Macierz nie spełnia założeń")
        return
    end

    L = eye(n); 
    D = zeros(n, 1); 

    for k = 1:n
        suma_dk = 0;
        for p = 1:k-1
            suma_dk = suma_dk + abs(L(k, p))^2 * D(p);
        end
        D(k) = A(k, k) - suma_dk;
        for i = k+1:n
            suma_lik = 0;
            for p = 1:k-1
                suma_lik = suma_lik + L(i, p) * conj(L(k, p)) * D(p);
            end
            L(i, k) = (A(i, k) - suma_lik) / D(k);
        end
    end
    D= diag(D);
end