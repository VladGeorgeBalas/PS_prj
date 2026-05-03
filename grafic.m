% Balas Vlad-George, 333AB
%--------------------------
% FILE: grafic.m
%
% FUNCTION: grafic
%
% CALL: grafic(B, A, xlines, ylines, titlu)
%
% Functie care deseneaza graficele caracteristicii unui filtru trimis prin numitor si numarator. 
% Plotarea fazei foloseste toate valorile mai mari decat eps din caracteristica in frecventa a
% filtrului. Rezolutia de plotare este 5000 conform instructiunilor.
%
% Argumente
%   B       - numitorul functiei de transfer
%   A       - numaratorul functiei de transfer
%   xlines  - [xlines_spectru; xlines_faza], in pulsatie normalizata
%   ylines  - [ylines_spectru; ylines_faza], va fi convertit in dB automat
%   titlu   - string care va fi afisat ca titlu al graficului.  
%
% Autor: Balas Vlad-George
% Creat: Decembrie, 2025

function [] = grafic(B, A, xlines, ylines, titlu)
rezolutie = 5000;

[H, w] = freqz(B, A, rezolutie);

figure_grafice = figure;
grafic_frecventa = subplot(2, 1, 1, 'Parent', figure_grafice);
grafic_faza = subplot(2, 1, 2, 'Parent', figure_grafice);
sgtitle(titlu, 'interpreter', 'latex');

axes(grafic_frecventa);
H_abs_db = db(abs(H));
plot(w, H_abs_db);
if ~isempty(xlines) && size(xlines,1) >= 1
    for x = xlines(1, :)
        xline(x * pi);
    end
end
if ~isempty(ylines) && size(ylines,1) >= 1
    for y = ylines(1, :)
        y_db = db(y);
        yline(y_db);
    end
end
ylabel("$|G(e^{\omega j})|\ (\ dB\ )$", 'interpreter', 'latex');
xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
xlim([0, pi]);
ylim([db(eps), 50]);                % limitez pana |H| < precizia numeria, acolo da erori numerice

axes(grafic_faza);
elemente_corecte = (abs(H) > eps);  % ca mai sus, doar ca pastrez doar valorile unde angle
                                    % lucreaza cu un abs(H) corect numeric
faza = unwrap(angle(H));
plot(w(elemente_corecte), faza(elemente_corecte));
if ~isempty(xlines) && size(xlines,1) >= 2
    for x = xlines(2, :)
        xline(x * pi);
    end
end
if ~isempty(ylines) && size(ylines,1) >= 2
    for y = ylines(2, :)
        y_db = db(y);
        yline(y_db);
    end
end
ylabel("$Faza\ (\ radian\ )$", 'interpreter', 'latex');
xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
xlim([0, pi]);

end