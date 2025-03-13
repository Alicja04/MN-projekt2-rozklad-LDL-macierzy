function drawplot(cholerr,errors, errorsSolution, condofA, n_min, n_max, param)
% drawplot Rysuje wykresy błędów i współczynnika uwarunkowania macierzy.
%   Funkcja rysuje wykresy błędów rozkładu, rozwiązania metody zaimplementowanej i metody bezpośredniej oraz
%   współczynnika uwarunkowania macierzy, z dwiema osiami Y. Parametr 
%   `param` decyduje, czy tytułem osi X będą rozmiary macierzy, czy kolejne macierze.
%
%   cholerr        - Błąd rozkładu
%   errors         - Błąd rozwiązania
%   errorsSolution - Błąd metody B\A
%   condofA        - Współczynnik uwarunkowania
%   n_min, n_max   - Zakres rozmiarów macierzy
%   param          - Typ osi X (2: rozmiary macierzy, inna wartość: kolejne macierze)

figure;

yyaxis left;
plot(n_min:n_max, cholerr, '-o', 'Color', 'b', 'LineWidth', 1.5); 
hold on;
plot(n_min:n_max, errors, '-o', 'Color', 'g', 'LineWidth', 1.5); 
plot(n_min:n_max, errorsSolution, '--o','Color', 'k', 'LineWidth', 1.5);
hold off;
ylabel('Błąd (err)');

yyaxis right;
plot(n_min:n_max, condofA, '-s', 'Color', 'r', 'LineWidth', 1.5);
ylabel('Cond macierzy');
ylim([0, max(condofA) * 1.1]);

if param == 2
    xlabel('Rozmiar macierzy');
    title('Błędy i współczynnik uwarunkowania macierzy w zależności od rozmiaru macierzy');
else
    xlabel('Kolejne macierze');
    title('Błędy i współczynnik uwarunkowania macierzy');
end
grid on;
legend({'Błąd rozkładu optymalnego', 'Błąd rozwiązania optymalnego','Błąd bezpośredniego rozw.', 'Cond'}, ...
       'Location', 'northwest');
end