% Balas Vlad-George, 333AB
%--------------------------
% FILE: criteriu.m
%
% FUNCTION: citeriu
%
% CALL: [scor_cost_mediu, scor_calitate_mediu, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, Delta_s)
% 
% Aplica criteriul de performanta cost/calitate pentru un filtru si tolerantele date. Returneaza scorul de cost, calitate
% cat si scorul compus obtinut cu cele doua.
%
% Argumente
%   b       - Numaratorul filtrului ce va fi evaluat
%   a       - Numaratorul filtrului ce va fi evaluat
%   Delta_p - Toleranta in banda de trecere
%   Delta_s - Toleranta in banda de stopare
%   omega_p - frecventa banda de trecere, in [0, pi]
%   omega_s - frecventa banda de stopare, in [0, pi]
%
% Returneaza
%   scor_cost_mediu     - Media aritmetica dintre scorul de cost in banda de trecere si de stopare
%   scor_calitate_mediu - Media aritmetica dintre scorul de calitate in banda de trecere si de stopare
%   scor_compus         - Scorul compus obtinut cu distanta euclidiana dintre scorul de cost mediu si cel de calitate mediu
%
% Foloseste:
%   TOL
%
% Autor: Balas Vlad-George
% Creat: Ianuarie, 2026

function [scor_cost_mediu, scor_calitate_mediu, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, Delta_s)
[H, w] = freqz(b, a, 5000);

% am facut maxim pentru a putea acomoda atat filtre IIR cat si filtre FIR
ordin = max(length(a), length(b)) - 1;

pass_band = w <= omega_p;
stop_band = omega_s <= w;

% in cazul in care filtrul nu se incadreaza in tolerante, primeste scorul 0 deoarece nu poate fi evaluat corect
% si nici nu respecta cerintele de proiectare
if ~tol(b, a, Delta_p, Delta_s, omega_p, omega_s)
    scor_cost_mediu = 0;
    scor_calitate_mediu = 0;
    scor_compus = 0;

    return;
end

% Criteriu de cost: cata calitate per pol avem
scor_cost_pass = (Delta_p - sum(abs(abs(H(pass_band)) - 1)) / nnz(pass_band)) / ordin / Delta_p;
scor_cost_stop = (Delta_s - sum(abs(abs(H(stop_band)))) / nnz(stop_band)) / ordin / Delta_s;
scor_cost_mediu = (scor_cost_pass + scor_cost_stop) / 2;

% Criteriu de calitate: cat de calitativ este rezultatul pentru ce tolerante avem
scor_calitate_pass = (Delta_p - max(abs(abs(H(pass_band)) - 1))) / Delta_p;
scor_calitate_stop = (Delta_s - max(abs(abs(H(stop_band))))) / Delta_s;
scor_calitate_mediu = 1/2 * (scor_calitate_stop + scor_calitate_pass);

scor_compus = sqrt(scor_cost_mediu ^ 2 + scor_calitate_mediu ^ 2);
end
