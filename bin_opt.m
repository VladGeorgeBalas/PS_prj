% Balas Vlad-George, 333AB
%--------------------------
% FILE: bin_opt.m
%
% FUNCTION: bin_opt
%
% CALL: [M, b, a] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M, metoda)
% 
% Gaseste cel mai mic ordin pentru ca filtrul FTJ sa fie in tolrante prin metoda indicata in indicatii.txt
% In cazul in care filtrul gasit este deja incadrat in tolerante, scade ordinul acestuia pana gaseste
% ordinul minim. In cazul in care filtrul nu este in tolerantele data, creste ordinul pana aceste se incadreaza.
%
% Argumente
%   Delta_p - Toleranta in banda de trecere
%   Delta_s - Toleranta in banda de stopare
%   omega_p - frecventa banda de trecere, in [0, pi]
%   omega_s - frecventa banda de stopare, in [0, pi]
%   Ts      - Perioada de esantionare
%   M       - Ordinul de la care pleaca gasirea ordinului minim
%   metoda  - metoda folosita pentru a gasi functia filtrului: 'ellip', 'firls', 'fir1', 'cheby1', 'cheby2'
%
% Returneaza
%   M       - Ordinul gasit
%   b       - numaratorul functiei de transfer gasita cu metoda data si M-ul returnat
%   a       - numitorul functiei de transfer gasita cu metoda data si M-ul returnat
%
% Foloseste:
%   TOL
%
% Autor: Balas Vlad-George
% Creat: Ianuarie, 2026

function [M, b, a] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M, metoda)

tustin = @(omega, t_s)((2 / t_s) * tan(omega / 2));
inverse_tustin = @(ohm, t_s)(2 * atan((t_s / 2) * ohm));

M_dat = M;

    function [b, a] = proiectare(Delta_p, Delta_s, omega_p, omega_s, Ts, M, M_dat, metoda)
        if strcmp(metoda, "fir1")
            M_p = 1 - Delta_p;

            omega_c = inverse_tustin(...
                tustin(omega_p, Ts) / nthroot((1 - M_p * M_p)/(M_p * M_p), 2 * M_dat),...
                Ts); % omega_c in radian

            b = fir1(M, omega_c / pi, "low");
            a = 1;
        elseif strcmp(metoda, "ellip")
            [b, a] = ellip(M, -20 * log10(1-Delta_p), -20 * log10(Delta_s), omega_p / pi);
        elseif strcmp(metoda, "firls")
            b = firls(M, [0, omega_p / pi, omega_s / pi, 1], [1, 1, 0, 0]);
            a = 1;
        elseif strcmp(metoda, "cheby1")
            [b, a] = cheby1(M, -20 * log10(1-Delta_p), omega_p / pi);
        elseif strcmp(metoda, "cheby2")
            [b, a] = cheby2(M, -20 * log10(Delta_s), omega_s / pi);
        end
    end

[b, a] = proiectare(Delta_p, Delta_s, omega_p, omega_s, Ts, M, M_dat, metoda);
if tol(b, a, Delta_p, Delta_s, omega_p, omega_s) == true
    while M > 0
        [b, a] = proiectare(Delta_p, Delta_s, omega_p, omega_s, Ts, M, M_dat, metoda);
        % tol(b, a, Delta_p, Delta_s, omega_p, omega_s)
        if tol(b, a, Delta_p, Delta_s, omega_p, omega_s) == true
            M = M - 1;
        else
            break;
        end
    end

    M = M + 1;
else
    while true
        [b, a] = proiectare(Delta_p, Delta_s, omega_p, omega_s, Ts, M, M_dat, metoda);
        % tol(b, a, Delta_p, Delta_s, omega_p, omega_s)
        if tol(b, a, Delta_p, Delta_s, omega_p, omega_s) == false
            M = M + 1;
        else
            break;
        end
    end
end

[b, a] = proiectare(Delta_p, Delta_s, omega_p, omega_s, Ts, M, M_dat, metoda);
end