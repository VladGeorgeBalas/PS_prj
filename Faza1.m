% Balas Vlad-George, 333AB

% [omega_p,omega_s,Delta_p,Ts] = PS_PRJ_3_Faza_1a(1,2)

% flag global de debug, cand rulez cu el true, printeaza niste
% valori de verificare, folosite pentru a gasi greseli cand ma verific
global DEBUG;
DEBUG = false;

% flag global de extra. cand este true, se mai scriu niste valori
% care ajuta cu orientarea si vederea pasilor in timpul functionarii.
global EXTRA;
EXTRA = false;

% Salvate dupa rularea functiei
omega_p = 0.9595
omega_s = 1.2633
Delta_p = 0.0623
Delta_s = Delta_p % in conformitate cu indicatiie date
Ts = 1.7027

% rezolutie
rez = 5000;

%% Functii utilitare
tustin = @(omega, t_s)((2 / t_s) * tan(omega / 2));
inverse_tustin = @(ohm, t_s)( 2 * atan((t_s / 2) * ohm));

pseudo_tustin = @(omega, t_s)((1 / t_s) * tan(omega / 2));
inverse_pseudo_tustin = @(ohm, t_s)(2 * atan(t_s * ohm));

%% Subpunctul a - Tustin

[B_tustin, A_tustin] = But_FTI(omega_p/ pi, omega_s / pi, Delta_p, Delta_s, Ts);
M = size(A_tustin, 2) - 1;
M_p = 1 - Delta_p;

% Calculare pulsatie de taiere
omega_c = inverse_tustin(...
    tustin(omega_p, Ts) / nthroot((1 - M_p * M_p)/(M_p * M_p), 2 * M),...
    Ts) % omega_c in radiani

% Matrice linii verticale
xlines = [omega_p / pi, omega_s / pi , omega_c / pi;
    omega_p / pi, omega_s / pi, omega_c / pi];

% Matrice linii orizontale
ylines = [1 + Delta_p, 1 - Delta_p, Delta_s;];

grafic(B_tustin, A_tustin, xlines, ylines, "Caracteristica in frecventa a filtrului cu $T_s = " + string(Ts) + "$ si $M = " + string(M) + "$");

% Debug, este aici doar sa compar ce am facut eu cu o functie oficiala
if DEBUG
    figure;
    freqz(B_tustin, A_tustin, 5000);
end

%% Subpunctul b - Pseudo Tustin

[B_pseudo_tustin, A_pseudo_tustin] = fnc(omega_p/ pi, omega_s / pi, Delta_p, Delta_s, Ts, "pseudo_tustin");
M = size(A_pseudo_tustin, 2) - 1;
M_p = 1 - Delta_p;

omega_c = inverse_pseudo_tustin(...
    pseudo_tustin(omega_p, Ts) / nthroot((1 - M_p * M_p)/(M_p * M_p), 2 * M),...
    Ts); % omega_c in radiani

% matrice linii verticale
xlines = [omega_p / pi, omega_s / pi , omega_c / pi;
    omega_p / pi, omega_s / pi, omega_c / pi];

% matrice linii orizontale
ylines = [1 + Delta_p, 1 - Delta_p, Delta_s;];

grafic(B_pseudo_tustin, A_pseudo_tustin, xlines, ylines, "Caracteristica in frecventa a filtrului cu $T_s = " + string(Ts) + "$ si $M = " + string(M) + "$");

% debug, este aici doar sa compar ce am facut eu cu o functie oficiala
if DEBUG
    figure;
    freqz(B_pseudo_tustin, A_pseudo_tustin, 5000);
end

[H_pseudo_tustin, w] = freqz(B_pseudo_tustin, A_pseudo_tustin, rez);
faza_pseudo_tustin = unwrap(angle(H_pseudo_tustin));

[H_tustin, w] = freqz(B_tustin, A_tustin, rez);
faza_tustin = unwrap(angle(H_tustin));

figure;
tiledlayout(2, 1);

nexttile;
H_err = abs(abs(H_tustin) - abs(H_pseudo_tustin));
norma_freventa = norm(H_err);

plot(w, H_err);
title("Norma este egala cu " + string(norma_freventa));

nexttile;
valori_corecte = ((abs(H_pseudo_tustin) > eps) & (abs(H_tustin) > eps));
faza_err = abs(faza_tustin(valori_corecte) - faza_pseudo_tustin(valori_corecte));
norma_faza = norm(faza_err);

plot(w(valori_corecte), faza_err);
title("Norma este egala cu " + string(norma_faza));


%% Subpunctul c

function grafic_c(t_s, B_cmp, A_cmp, rez, omega_p, omega_s, Delta_p, Delta_s)
% Descriere
%   Mica functie pe care o refolosesc sa desenez cele 2 matrici de grafice
%
% Argumente
%   t_s     - O matrice de forma [t_s_1, t_s_2, t_s_3, t_s_4] care contine cele 4 perioade de esantionare
%   B_cmp   - Numitorul functiei de transfer cu care comparam
%   A_cmp   - Numaratorul functiei de transfer cu care comparam
%   rez     - Rezolutia pentru freqz

%% Functii utilitare

tustin = @(omega, t_s)((2 / t_s) * tan(omega / 2));
inverse_tustin = @(ohm, t_s)( 2 * atan((t_s / 2) * ohm));

%% Flag-uri globale

global EXTRA;
global DEBUG;

%% Corpul functiei

figure;

[H_cmp, w] = freqz(B_cmp, A_cmp, rez);
faza_cmp = unwrap(angle(H_cmp));

for i = 1:4
    [B, A] = fnc(omega_p/ pi, omega_s / pi, Delta_p, Delta_s, t_s(i), "tustin");
    [H, w] = freqz(B, A, rez);

    if EXTRA
        norm(A - A_cmp)
        norm(B - B_cmp)
    end

    M = size(A, 2) - 1;
    M_p = 1 - Delta_p;

    omega_c = inverse_tustin(...
        tustin(omega_p, t_s(i)) / nthroot((1 - M_p * M_p)/(M_p * M_p), 2 * M),...
        t_s(i)); % omega_c in radiani

    valori_corecte = (abs(H_cmp) > eps & abs(H) > eps);

    %% Plot caracteristica frecventa
    subplot(4, 4, i);

    H_abs_db = db(abs(H));
    plot(w, H_abs_db);

    linii_verticale = [omega_p, omega_c, omega_s];
    linii_orizontale = [1 + Delta_p, 1 - Delta_p, Delta_s];
    linii_orizontale = db(linii_orizontale);

    xline(linii_verticale);
    yline(linii_orizontale);

    ylabel("$|G(e^{\omega j})|\ (\ dB\ )$", 'interpreter', 'latex');
    xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
    xlim([0, pi]);
    ylim([db(eps), 50]);
    title(sprintf('$T_{s}=$%.3e', t_s(i)), 'Interpreter','latex');

    %% Plot grafic eroare
    subplot(4, 4, 4 + i);

    H_err = abs(abs(H_cmp) - abs(H));
    norma_freventa = norm(H_err);

    plot(w, H_err);
    title(sprintf('Norma = %.3e', norma_freventa), 'Interpreter','latex');

    %% Plot caracteristica in faza

    subplot(4, 4, 8 + i);
    faza = unwrap(angle(H));
    plot(w(valori_corecte), faza(valori_corecte));

    xline(linii_verticale);

    ylabel("$Faza\ (\ radian\ )$", 'interpreter', 'latex');
    xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
    xlim([0, pi]);

    %% Plot grafic eroare

    subplot(4, 4, 12 + i);
    faza_err = abs(faza_cmp(valori_corecte) - faza(valori_corecte));
    norma_faza = norm(faza_err);

    plot(w(valori_corecte), faza_err);
    title(sprintf('Norma = %.3e', norma_faza), 'Interpreter','latex');
end
end

grafic_c([0.1 * Ts, Ts / 4, Ts / 2, 3 * Ts / 4], B_tustin, A_tustin, rez, omega_p, omega_s, Delta_p, Delta_s);
grafic_c([5 * Ts / 4, 7 * Ts / 4, 9 * Ts / 4, 3 * Ts], B_tustin, A_tustin, rez, omega_p, omega_s, Delta_p, Delta_s);

%% Subpuntul d

function R = grafic_d(Ts, rez, omega_p, omega_s, delta_combo)
% Descriere
%   Mica functie pe care o refolosesc sa desenez cele 2 matrici de grafice
%
% Argumente
%   delta_combo - [[delta_p, delta_s]; [delta_p, delta_s]; etc]
%   rez     - Rezolutia pentru freqz

%% Functii utilitare

tustin = @(omega, t_s)((2 / t_s) * tan(omega / 2));
inverse_tustin = @(ohm, t_s)( 2 * atan((t_s / 2) * ohm));

%% Flag-uri globale

global EXTRA;
global DEBUG;

%% Corpul functiei

R = [];
figure;

for i = 1:size(delta_combo, 1)
    [B, A] = fnc(omega_p/ pi, omega_s / pi, delta_combo(i, 1), delta_combo(i, 2), Ts, "tustin");
    [H, w] = freqz(B, A, rez);

    M = size(A, 2) - 1;
    M_p = 1 - delta_combo(i, 1);

    R = [R; delta_combo(i, :), M];

    omega_c = inverse_tustin(...
        tustin(omega_p, Ts) / nthroot((1 - M_p * M_p)/(M_p * M_p), 2 * M),...
        Ts); % omega_c in radiani

    valori_corecte = (abs(H) > eps);

    %% Plot caracteristica frecventa
    subplot(4, 4, floor(i / 5) * 8 + mod(i - 1, 4) + 1);

    H_abs_db = db(abs(H));
    plot(w, H_abs_db);

    linii_verticale = [omega_p, omega_c, omega_s];
    linii_orizontale = [1 + delta_combo(i, 1), 1 - delta_combo(i, 1), delta_combo(i, 2)];
    linii_orizontale = db(linii_orizontale);

    xline(linii_verticale);
    yline(linii_orizontale);

    ylabel("$|G(e^{\omega j})|\ (\ dB\ )$", 'interpreter', 'latex');
    xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
    xlim([0, pi]);
    ylim([db(eps), 50]);

    title({ ...
        sprintf('M = %d', M), ...
        sprintf('$\\Delta_p = %.4g\\quad \\Delta_s = %.4g$', ...
        delta_combo(i,1), delta_combo(i,2)) ...
        }, 'Interpreter','latex');

    %% Plot caracteristica in faza

    subplot(4, 4, floor(i / 5) * 8 + mod(i - 1, 4) + 1 + 4);
    faza = unwrap(angle(H));
    plot(w(valori_corecte), faza(valori_corecte));

    xline(linii_verticale);

    ylabel("$Faza\ (\ radian\ )$", 'interpreter', 'latex');
    xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
    xlim([0, pi]);
end
end

delta_combo = [
    [Delta_p/ 2, Delta_p / 2];
    [Delta_p/ 2, Delta_p];
    [Delta_p/ 2, 3 * Delta_p / 2];
    [Delta_p/ 2, 2 * Delta_p];

    [Delta_p, Delta_p / 2];
    [Delta_p, Delta_p];
    [Delta_p, 3 * Delta_p / 2];
    [Delta_p, 2 * Delta_p];
    ];
R1 = grafic_d(Ts, rez, omega_p, omega_s, delta_combo);

delta_combo = [
    [3 * Delta_p/ 2, Delta_p / 2];
    [3 * Delta_p/ 2, Delta_p];
    [3 * Delta_p/ 2, 3 * Delta_p / 2];
    [3 * Delta_p/ 2, 2 * Delta_p];

    [2 * Delta_p, Delta_p / 2];
    [2 * Delta_p, Delta_p];
    [2 * Delta_p, 3 * Delta_p / 2];
    [2 * Delta_p, 2 * Delta_p];
    ];
R2 = grafic_d(Ts, rez, omega_p, omega_s, delta_combo);

if EXTRA
    % un grafic care arata cum depinde M de delta_p si delta S
    % este sub flag-ul EXTRA, deoarece nu este cerut in cerinta
    % dar arata foarte clar cum scade, si mi s-a parut interesant
    figure;
    scatter3([R1(:, 1); R2(:, 1)], [R1(:, 2); R2(:, 2)], [R1(:, 3); R2(:, 3)], 60, [R1(:, 3); R2(:, 3)], 'filled');
    grid on;

    xlabel('\Delta_p');
    ylabel('\Delta_s');
    zlabel('M (ordin)');
    title('Scăderea ordinului M când cresc toleranțele \Delta_p, \Delta_s');
end

%% Subpunctul e

function grafic_e(b_1, a_1, b_2, a_2, omega_p, omega_s, omega_c, Delta_p, Delta_s, rez, M, tip_filtru)
% Descriere
%   Functie care ploteaza in paralele graficul raspunsului in fraventa/faza al functie de transfer 1(b_1 si a_1) si eroarea fata de functia de transfer 2
%
% Argumente
%   rez     - Rezolutia pentru freqz

%% Functii utilitare

%% Flag-uri globale

global EXTRA;
global DEBUG;

%% Corpul functiei

figure;

[H_cmp, w] = freqz(b_2, a_2, rez);
faza_cmp = unwrap(angle(H_cmp));

[H, w] = freqz(b_1, a_1, rez);
faza = unwrap(angle(H));

valori_corecte = (abs(H_cmp) > eps & abs(H) > eps);

%% Plot caracteristica frecventa
subplot(2, 2, 1);

H_abs=abs(H);
plot(w, H_abs);

linii_verticale = [omega_p, omega_c, omega_s];
linii_orizontale = [1 + Delta_p, 1 - Delta_p, Delta_s];

xline(linii_verticale);
yline(linii_orizontale);

ylabel("$|G(e^{\omega j})|$", 'interpreter', 'latex');
xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
xlim([0, pi]);
ylim([0, 1.5]);
title(sprintf(tip_filtru + " cu ordinul "+ '$M = %d$', M), 'Interpreter','latex');

%% Plot grafic eroare
subplot(2, 2, 2);

H_err = abs(abs(H_cmp) - abs(H));
norma_freventa = norm(H_err);

plot(w, H_err);
title(sprintf('Norma = %.3e', norma_freventa), 'Interpreter','latex');

%% Plot caracteristica in faza

subplot(2, 2, 3);
plot(w(valori_corecte), faza(valori_corecte));

xline(linii_verticale);

ylabel("$Faza\ (\ radian\ )$", 'interpreter', 'latex');
xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
xlim([0, pi]);

%% Plot grafic eroare

subplot(2, 2, 4);
faza_err = abs(faza_cmp(valori_corecte) - faza(valori_corecte));
norma_faza = norm(faza_err);

plot(w(valori_corecte), faza_err);
title(sprintf('Norma = %.3e', norma_faza), 'Interpreter','latex');
end

% Calculare pulsatie de taiere
if DEBUG
    disp(sprintf('M = d%', M))
end
omega_c = inverse_tustin(...
    tustin(omega_p, Ts) / nthroot((1 - M_p * M_p)/(M_p * M_p), 2 * M),...
    Ts); % omega_c in radian

[M_fir1, ~, ~] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M, "fir1");
% nu scadem din M_fir1 un 1 deoarece functia bin_opt returneaza ordinul real si ultimul corect
b = fir1(M_fir1, omega_c / pi, "low");
grafic_e(b, 1, B_tustin, A_tustin, omega_p, omega_s, omega_c, Delta_p, Delta_s, rez, M_fir1, "FIR obtinut cu metoda ferestrei")


[M_firls, ~, ~] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M, "firls");
% nu scadem din M_firls un 1 deoarece functia bin_opt returneaza ordinul real si ultimul corect
b = firls(M_firls, [0, omega_p / pi, omega_s / pi, 1], [1, 1, 0, 0]);
grafic_e(b, 1, B_tustin, A_tustin, omega_p, omega_s, omega_c, Delta_p, Delta_s, rez, M_firls, "FIR obtinut cu metoda celor mai mici patrate")


