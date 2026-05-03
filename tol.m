% Balas Vlad-George, 333AB
%--------------------------
% FILE: tol.m 
% 
% FUNCTION: tol
%
% CALL: ok = tol(b, a, Delta_p, Delta_s, omega_p, omega_s)
% 
% Functia returneaza o true daca filtrul se incadreaza in tolerante la o rezolutie de 5000.
% In cazul in care nu se incadreaza, funtia returneaza false.
%
% Argumente
%   b       - Numaratorul filtrului ce va fi evaluat
%   a       - Numaratorul filtrului ce va fi evaluat
%   Delta_p - Toleranta in banda de trecere
%   Delta_s - Toleranta in banda de stopare
%   omega_p - frecventa banda de trecere, in [0, pi]
%   omega_s - frecventa banda de stopare, in [0, pi]
%
% Iesire
%   ok      - true daca filtrul se incadreaza, false daca nu
%
% Autor: Balas Vlad-George
% Creat: Ianuarie, 2026

function ok = tol(b, a, Delta_p, Delta_s, omega_p, omega_s)
    rez = 5000; % rezolutia 5000, conform indicatiilor
    [H, w] = freqz(b, a, rez);

    omega_p_disc = find(w <= omega_p, 1, 'last');
    omega_s_disc = find(w >= omega_s, 1, 'first');

    H_abs_p = abs(H(omega_p_disc));
    H_abs_s = abs(H(omega_s_disc));

    ok = (H_abs_p >= (1 - Delta_p)) && (H_abs_p <= (1 + Delta_p)) && (H_abs_s <= Delta_s);
end
