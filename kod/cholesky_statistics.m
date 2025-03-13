function resultTable = cholesky_statistics(A, B, Z)
% cholesky_statistics Oblicza wskaźniki statystyczne dla równania XA = B
% przy użyciu metod rozkładu Cholesky'ego.
%
% Wejście:
%   A - Macierz współczynników
%   B - Macierz prawej strony
%   Z - Wektor oczekiwanych rozwiązań
%
% Wyjście:
%   resultTable - Tabela z wynikami (błąd, stabilność, poprawność, błąd rozkładu)
%
% Metody:
%   1. rozklad_cholesky_LDLh_optimized
%   2. matrix_decomposition
%   3. ldl
%   4. B/A

  if ~ifProperMatrix(A) 
    disp("Macierz nie spełnia założeń")
    return
  end
  if size(B,2) ~= size(A,2)
    disp("Niezgodne wymiary")
    return
  end

  disp("Współczynnik uwarunkowania macierzy A");
  disp(cond(A));

  names = {'błąd względny rozwiązania', 'współczynnik stabilności', 'współczynnik poprawności', 'błąd rozkładu'};
  
  resultData = zeros(4, 4);  % 4 metody, 4 wskaźniki
  
  %metoda 2 (rozklad_cholesky_LDLh_optimized)
  [L2, D2] = rozklad_cholesky_LDLh_optimized(A);
  AChol2 = L2 * D2 * (L2');
  err2 = norm(A - AChol2) / norm(A);
  x = solveUsingCholesky(A, B, 2);
  resultData(1, 1) = norm(x - Z) / norm(Z);
  resultData(1, 2) = resultData(1, 1) / cond(A);
  resultData(1, 3) = norm(B - x * A) / (norm(A) * norm(x));
  resultData(1, 4) = err2;  % Błąd rozkładu dla metody 2

  %metoda 5 (matrix_decomposition)
  [L5, D5] = matrix_decomposition(A);
  AChol5 = L5 * D5 * (L5');
  err5 = norm(A - AChol5) / norm(A);
  x = solveUsingCholesky(A, B, 4);
  resultData(2, 1) = norm(x - Z) / norm(Z);
  resultData(2, 2) = resultData(2, 1) / cond(A);
  resultData(2, 3) = norm(B - x * A) / (norm(A) * norm(x));
  resultData(2, 4) = err5;  % Błąd rozkładu dla metody 5

  %metoda 3 (ldl - wbudowana funkcja)
  [L3, D3] = ldl(A);
  AChol3 = L3 * D3 * (L3');
  err3 = norm(A - AChol3) / norm(A);
  x = solveUsingCholesky(A, B, 3);
  resultData(3, 1) = norm(x - Z) / norm(Z);
  resultData(3, 2) = resultData(3, 1) / cond(A);
  resultData(3, 3) = norm(B - x * A) / (norm(A) * norm(x));
  resultData(3, 4) = err3;  % Błąd rozkładu dla metody 3

  %metoda 4 (bezpośrednie rozwiązanie B/A)
  x = B / A;
  resultData(4, 1) = norm(x - Z) / norm(Z);
  resultData(4, 2) = resultData(4, 1) / cond(A);
  resultData(4, 3) = norm(B - x * A) / (norm(A) * norm(x));
  resultData(4, 4) = -1;  % Dla metody B/A, błąd rozkładu = -1

  % Tworzenie tabeli wyników
  resultTable = array2table(resultData, 'VariableNames', names, 'RowNames', ...
                            {'rozklad_cholesky_LDLh_optimized', 'matrix_decomposition', 'ldl (wbudowane)', 'B/A(bezpośrednie rozwiązanie)'});
  
  % Wyświetlenie wyników
  disp("Tabela wyników:");
  disp(resultTable);
end
