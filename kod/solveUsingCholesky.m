function [wynik, L, D] = solveUsingCholesky(A, B, param)
%solveUsingCholesky Funkcja rozwiązuje równanie macierzowe postaci XA = B
%   
% Dane wejściowe:
%   * A - macierz rozkładana w celu rozwiązania równania (spełniająca warunki : hermitowska, pięciodiagonalna, dodatnio określona)
%   * B - macierz , której kolejne wiersze to kolejne wektory b równania XA=b
%
% Dane wyjściowe:
%   * wynik - macierz, której kolejne wiersze to wyniki dla kolejnych
%     wierszy macierzy B

    if size(B,2) ~= size(A,2)
        disp("Niezgodne wymiary")
        return
    end
    
    if ~ifProperMatrix(A) 
        disp("Macierz nie spełnia założeń")
        return
    end
    
    if param == 1
        [L,D] =rozklad_cholesky_LDLh(A); %niedokładna meotda, nieużyta w projekcie
    elseif param == 2    
        [L,D] = rozklad_cholesky_LDLh_optimized(A); %ok
    elseif param == 3
        [L,D] = ldl(A);% wbudowana
    elseif param == 4
        [L,D] = matrix_decomposition(A); %ok
    else
        wynik = B/A; % wbudowana
        return 
    end
    
        Z = B / L';
        Y = Z / D;
        X = Y / L;
        wynik = X;
   
end