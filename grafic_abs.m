function [] = grafic_abs(B, A, xlines, ylines, titlu)
% GRAFIC - functie care deseneaza graficele caracteristicii unui filtru in
%   functie de numitor si numarator
%
% Argumente
%   B   -   numitorul functiei de transfer
%   A   -   numaratorul functiei de transfer
%   xlines  -   [xlines_spectru; xlines_faza], in pulsatie normalizata
%   ylines  -   [ylines_spectru; ylines_faza], va fi convertit in dB automat

% rezolutie 5000, conform indicatiilor 
rezolutie = 5000;

[H, w] = freqz(B, A, rezolutie);

figure_grafice = figure;
grafic_frecventa = subplot(2, 1, 1, 'Parent', figure_grafice);
grafic_faza = subplot(2, 1, 2, 'Parent', figure_grafice);
sgtitle(titlu, 'interpreter', 'latex');

axes(grafic_frecventa);
H_abs = abs(H);
plot(w, H_abs);
if ~isempty(xlines) && size(xlines,1) >= 1
    for x = xlines(1, :)
        xline(x * pi);
    end
end
if ~isempty(ylines) && size(ylines,1) >= 1
    for y = ylines(1, :)
        yline(y);
    end
end
ylabel("$|G(e^{\omega j})|$", 'interpreter', 'latex');
xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
xlim([0, pi]);
ylim([0, 1.5]);                % limitez pana H < precizia numeria, acolo da erori numerice

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
        yline(y);
    end
end
ylabel("$Faza\ (\ radian\ )$", 'interpreter', 'latex');
xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
xlim([0, pi]);

end