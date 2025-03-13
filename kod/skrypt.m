%% Dostępne Z do wyboru
Z = 100* randn(100, 5) + 100i * randn(100, 5);
Z2 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 1 1 1 1;
    1 2 3 4 5; 6 7 8 9 10; 100 -100 200 -300 20;
    -1 1 -1 1 -1; 0 1 -1 1 -1; 100 -100 100 -100 0;
    4 80 -4 80 100; 1i 1i 1i 1i 1i;
    2i 3i 4i 5i 6i; 0 1i -1i 1i -1i; 1i 1 1i 1 -1i];
%% Ustawienie Z
%miejsce na ustawienie Z używanego w większości przykładów poniżej (Z
%odpowiada za losowo generowany wektor, Z2 za podany wektor o widocznych
%wartościach)
    Z = Z;
%% Macierze rzeczywiste
A0 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
A1 =[5 3 0 0 0; 3 6 2 0 0; 0 2 7 1 0; 0 0 1 10 0; 0 0 0 0 1];
A2 = [4 -1 0.5 0 0;-1 4 -1 0.5 0;0.5 -1 4 -1 0.5; 0 0.5 -1 4 -1; 0 0 0.5 -1 4];
A3 = [4 -1 0 0 0;-1 4 -1 0 0;0 -1 4 -1 0; 0 0 -1 4 -1; 0 0 0 -1 4];
A4  =[[ 6.9378 -2.0657 -1.4554  0.0000  0.0000]
 [-2.0657  4.3818 -1.2078 -0.7302  0.0000]
 [-1.4554 -1.2078  8.1936 -2.5097 -0.4876]
 [ 0.0000 -0.7302 -2.5097  3.5736 -1.6042]
 [ 0.0000  0.0000 -0.4876 -1.6042  4.4132]];

A5 = [[ 23.2904  -9.4886  -5.8119   0.0000   0.0000]
 [ -9.4886  40.1439  -8.6689  -9.2368   0.0000]
 [ -5.8119  -8.6689  19.3262 -10.5067  -4.4213]
 [  0.0000  -9.2368 -10.5067  34.1298  -7.6254]
 [  0.0000   0.0000  -4.4213  -7.6254  10.6095]];
A6 = [[ 28.1750 -21.4240  -7.6822   0.0000   0.0000]
 [-21.4240  87.9700 -19.2470  -9.0065   0.0000]
 [ -7.6822 -19.2470  33.8959 -17.0160  -9.3408]
 [  0.0000  -9.0065 -17.0160  33.8335 -17.9227]
 [  0.0000   0.0000  -9.3408 -17.9227  68.6256]];
matrix = cell({A0; A1; A2; A3; A4; A5; A6});
errors = zeros(7, 1); 
condofA = zeros(7, 1);
cholerr = zeros(7, 1);
errorsSolution = zeros(7, 1);

for n = 1:7
    A = matrix{n};
    condofA(n) = cond(A);
    B = Z * A;  
    [x,L,D] = solveUsingCholesky(A, B, 2);  
    err = norm(x - Z) / norm(Z);
    errors(n) = err;
    cholerr(n) = norm(A-L*D*L')/norm(A);
    errorsSolution(n) = norm(solveUsingCholesky(A,B,5)-Z)/norm(Z);
end
drawplot(cholerr, errors, errorsSolution, condofA, 1, 7, 1);

%dla macierzy rzeczywistych , ale Z dowolne rozwiązania dokładne lub
%bliskie dokładnym błąd rzedu 10^(-15), dla Z wygenerowanego losowo wyniki
%trochę gorsze, ale dalej w zakresie tolerancji (3*10^(-15)), błąd rośnie
%ze wzrostem współczynnika uwarunkowania macierzy cond(A)
%jednostkowa dokłądne rozwiązania
%% Macierze zespolone
A7 = [10 1-1i 0 0 0;
		  1+1i 5 0 0 0; 
		  0 0 6 0 0; 
		  0 0 0 7 0; 
		  0 0 0 0 7];
A8 = A1 + [0 1i 0 0 0; -1i 0 1/2i 0 0; 0 -1/2i 0 1/3i 0; 0 0 -1/3i 0 1/4i; 0 0 0 -1/4i 0];
A9 = [2 -1i 0 0 0; 1i 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
A10 = [ 4 1+1i 0 0 0; 1-1i 3 0 0 0; 0 0 5 1+1i 0; 0 0 1-1i 6 0; 0 0 0 0 7]; % wszystkie powyżej to tylko trzy przekątne niezerowe
A11 = [4 1+1i 2 0 0; 1-1i 3 0 2+1i 0; 2 0 5 1+1i -1; 0 2-1i 1-1i 6 0; 0 0 -1 0 7];
A12 = A10 + [0 0 0 0 0; 0 0 1+1i 0 0; 0 1-1i 0 0 0; 0 0 0 0 1+1i; 0 0 0 1-1i 0];
A13 = A4+[0 1i 2i 0 0; -1i 0 1/2i 3i 0; -2i -1/2i 0 1/3i 4i; 0 -3i -1/3i 0 1/4i; 0 0 -4i -1/4i 0] + [2 0 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 2];

matrixComplex = cell({A7; A8; A9; A10; A11; A12; A13});
errorsComplex = zeros(7, 1);
condofAComplex = zeros(7, 1);
cholerrComplex = zeros(7, 1);
errorsSolutionComplex = zeros(7, 1);

for n = 1:7
    A = matrixComplex{n};
    condofAComplex(n) = cond(A);
    B = Z * A; 
    [x,L,D] = solveUsingCholesky(A, B, 2); 
    err = norm(x - Z) / norm(Z);
    errorsComplex(n) = err;
    cholerrComplex(n) = norm(A-L*D*L')/norm(A);
    errorsSolutionComplex(n) = norm(solveUsingCholesky(A,B,5)-Z)/norm(Z);
end
drawplot(cholerrComplex, errorsComplex, errorsSolutionComplex, condofAComplex, 1, 7, 1);

%dla macierzy o dużym cond nawet niewielkie błędy rozkładu (rzędu 10^(-16))
%mogą mieć wpływ na dokładność rozwiązania, rozwiązania tych układów też
%mieszczą się w błędzie tolerancji, ale są większe od bezpośredniego
%rozwiązania B\A (kótre jest dokładne), nawet dla dokładnego rozkładu bład
%rozwiązania jest rzedu 10^-16

%% Macierze zespolone o bardzo dużym współczynniku uwarunkowania
A14 = A6 + [0 10i 2i 0 0; -10i 0 1/2i 3i 0; -2i -1/2i 0 1/3i 4i; 0 -3i -1/3i 0 1/4i; 0 0 -4i -1/4i 0];
A15 = A5 + [0 1i 2i 0 0; -1i 0 1/20i 3i 0; -2i -1/20i 0 1/3i 4i; 0 -3i -1/3i 0 1/4i; 0 0 -4i -1/4i 0]; %macierz o bardzo dużym cond

cholesky_statistics(A14, Z*A14, Z);
cholesky_statistics(A15, Z*A15, Z);
%uzywając zoptymalizowanego algorytmu dostajemy trochę większe błędy
%rozwiązań (ale im większy cond tym ta różnica będzie bardziej zanikać)
%zyskujemy trochę lepsza stabilność, dlatego jeśli nie zależy nam na jak
%najbardziej dokładnym wyniku, to warto skorzystać z opcji optimized, ze
%wzgklędu na lepszą złożoność algorytmu
%% Macierze większych rozmiarów (6x6 i 8x8)
A16 = [ 4  -1   -1/2+1i  0   0    0    0    0; 
     -1   5  -1   -1/2   0   0    0    0; 
      -1/2-1i  -1   6  -1    -1/2   0   0    0; 
     0  -1/2  -1   7   -1    -1/2   0   0; 
      0  0  -1/2  -1    7   -1    -1/2   0; 
      0   0  0 -1/2   -1    6   -1    -1/2; 
      0   0   0  0   -1/2   -1    5   -1; 
      0   0   0   0   0   -1/2   -1    4];
Z8 = 100* randn(100, 8) + 100i * randn(100, 8);
cholesky_statistics(A16, Z8*A16, Z8);
A17 = [5 -2 1 0 0 0;
    -2 6 -2 1 0 0;
    1 -2 7 -2 1 0;
    0 1 -2 8 -2 1;
    0 0 1 -2 9 -2;
    0 0 0 1 -2 10];
Z6 = 100* randn(100, 6) + 100i * randn(100, 6);
cholesky_statistics(A17, Z6*A17, Z6);
% niewielkie błędy ale cond też są niewielkie około 4, różnice pomiędzy
% funkcjami nie są wielkie
%% Macierze kolejnych rozmiarów generowane losowo
n_max = 100;
errorsn = zeros(n_max-4, 1); % wektor do przechowywania błędów
condofAn = zeros(n_max-4, 1);
cholerrn = zeros(n_max-4, 1);
errorsSolutionn = zeros(n_max-4, 1);
matrixn = generate_positive_definite_pentadiagonal(n_max);

for n = 5:n_max
    A = matrixn{n - 4};
    condofAn(n-4) = cond(A);
    Z = 100* randn(100, n) + 100i * randn(100, n);
    B = Z * A;  
    [x,L,D] = solveUsingCholesky(A, B, 2); 
    err = norm(x - Z) / norm(Z);  
    errorsn(n - 4) = err;
    cholerrn(n-4) = norm(A-L*D*L')/norm(A);
    errorsSolutionn(n-4) = norm(solveUsingCholesky(A,B,5)-Z)/norm(Z);
end

drawplot(cholerrn, errorsn, errorsSolutionn, condofAn, 5, n_max, 2);