function [booleanValue] = ifProperMatrix(A)
%ifProperMatrix Funkcja sprawdza czy macierz wejściowa A
%               spełnia poniższe warunki : 
%               - jest dodatnio określona
%               - jest hermitowska
%               - jest pięciodiagonalna
%   Dane wejściowe:
%   - A - macierz kwadratowa, dla której sprawdzamy warunki
%
%   Dane wyjściowe:
%   - booleanValue - wartość true/false oznajmiająca czy macierz wejściowa
%   spełnia warunki
%   

    booleanValue = false;
    [n, m] = size(A);
    if n ~= m 
        disp("Macierz nie jest kwadratowa");
        return
    end
    if n < 5 
        disp("Macierz nie jest pięciodiagonalna");
        return
    end
   
    if ~isequal(A,A')
        disp("Macierz nie jest hermitowska");
        return
    end
  
    for i = 3:n
        if any(diag(A,i)~=-0) || any(diag(A,-i) ~= 0)
            disp("Macierz nie jest pięiodiagonalna (posiada niezerowe elementy poza wyznaczonymi przekątnymi)")
            return
        end
    end
   
    wartosciWlasne = eig(A);
    if all(wartosciWlasne > 0) %prawdziwe dla macierzy hermitowskich
        booleanValue = true;
    else
        disp("Macierz nie jest dodatnio określona (nie wszystkie wartości własne są dodatnie)");
    end
end