function matrices = generate_positive_definite_pentadiagonal(n)
    % Generuje pięciodiagonalne macierze dodatnio określone od 5x5 do nxn.
    % Wejście:
    %   n - maksymalny rozmiar macierzy (n >= 5)
    % Wyjście:
    %   matrices - komórka zawierająca macierze pięciodiagonalne

    if n < 5
        error('Rozmiar macierzy musi być co najmniej 5x5');
    end
    matrices = cell(n - 4, 1);
    for size = 5:n
        A = zeros(size);
        for i = 1:size
            A(i, i) = 1 + 100*rand(); 
            if i > 1
                A(i, i-1) = -1 + 0.1 * rand() + 10i*rand();
                A(i-1, i) = conj(A(i, i-1));
            end
            if i > 2
                A(i, i-2) = -0.5 + 0.1 * rand() + 10i*rand(); 
                A(i-2, i) = conj(A(i, i-2));
            end
        end
        A = A + size * eye(size);
        matrices{size - 4} = A;
    end
end
