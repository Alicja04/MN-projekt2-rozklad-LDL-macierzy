# Rozwiązywanie macierzowego układu równań metodą Cholesky’ego-Banachiewicza

## Opis projektu
Projekt realizowany w ramach przedmiotu **Metody numeryczne** na **Politechnice Warszawskiej** (wykonany w roku akademickim 2024/2025). Celem projektu jest analiza skuteczności i zastosowania **metody Cholesky’ego-Banachiewicza** w rozwiązywaniu macierzowego układu równań $X*A=B$ w dziedzinie zespolonej poprzez rozkład $LDL^h$ macierzy A.

**Język programowania**: MATLAB

**Autor**: Alicja Przeździecka 

## Kluczowe obserwacje
- **Bardzo duża dokładność** uzyskiwana jest dla macierzy rzeczywistych i zespolonych z maksymalnie trzema niezerowymi przekątnymi, co czyni metodę skuteczną dla rzadkich układów.
- **Lepsze uwarunkowanie macierzy** prowadzi do mniejszych błędów numerycznych, co poprawia stabilność wyników.
- **Zaimplementowana funkcja jest szczególnie przydatna dla dużych macierzy**, gdzie klasyczne metody eliminacyjne mogą wymagać zbyt dużej ilości pamięci i czasu obliczeniowego.
- **Metody oparte na rozkładzie $LDL^h$ zapewniają wyższą stabilność numeryczną**, ale mogą być wrażliwe na bardzo źle uwarunkowane macierze.
