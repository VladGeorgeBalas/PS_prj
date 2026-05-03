% Balas Vlad-George, 333AB

% omega_p = PS_PRJ_3_Faza_3(1, 2);

% Salvate dupa rularea functiei
omega_p = 1.4774
omega_s = omega_p + pi/33
Delta_p = 0.05 % 5% conform indicatiei
Ts = 1.7027

inversa_db = @(val_db)(10^(val_db/(20)));

% toate variabilele ce le ca sa plotez graficul cu criteriul de cost si calitate
ord_but = [];
ord_elip = [];
ord_cheby1 = [];
ord_cheby2 = [];

delta_s_but = [];
delta_s_elip = [];
delta_s_cheby1 = [];
delta_s_cheby2 = [];

scor_calitate_but = [];
scor_calitate_elip = [];
scor_calitate_cheby1 = [];
scor_calitate_cheby2 = [];

scor_cost_but = [];
scor_cost_elip = [];
scor_cost_cheby1 = [];
scor_cost_cheby2 = [];

scor_compus_but = [];
scor_compus_elip = [];
scor_compus_cheby1 = [];
scor_compus_cheby2 = [];

scor_max_but = [0, 0];
scor_max_elip = [0, 0];
scor_max_cheby1 = [0, 0];
scor_max_cheby2 = [0, 0];

ds = [];

for i = linspace(db2mag(-100), db2mag(-30), 1000)
    % aplicare criteriu pt filtrul Butterworth
    [b, a] = fnc(omega_p / pi, omega_s / pi, Delta_p, i, Ts, "tustin");
    n_but = size(a, 2) - 1;
    ord_but = [ord_but, n_but];

    [scor_cost, scor_calitate, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, i);
    scor_calitate_but = [scor_calitate_but, scor_calitate];
    scor_cost_but = [scor_cost_but, scor_cost];
    scor_compus_but = [scor_compus_but, scor_compus];

    delta_s_but = [delta_s_but, i];

    if scor_compus > scor_max_but(1)
        scor_max_but = [max(scor_max_but(1), scor_compus), i];
    end

    % aplicare criteriu pt filtrul eliptic cu ordin minim
    % folosesc funtia data de MATLAB pentru gasirea ordinului minim
    [n_elip, Wn] = ellipord(omega_p / pi, omega_s / pi, -20 * log10(1 - Delta_p), -20 * log10(i));
    if n_elip~=inf
        ord_elip = [ord_elip, n_elip];
        [b, a] = ellip(n_elip, -20 * log10(1 - Delta_p), -20 * log10(i), Wn);

        [scor_cost, scor_calitate, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, i);
        scor_cost_elip = [scor_cost_elip, scor_cost];
        scor_calitate_elip = [scor_calitate_elip, scor_calitate];
        scor_compus_elip = [scor_compus_elip, scor_compus];

        delta_s_elip = [delta_s_elip, i];

        if scor_compus > scor_max_elip(1)
            scor_max_elip = [max(scor_max_elip(1), scor_compus), i];
        end
    end

    % folosesc funtia data de MATLAB pentru gasirea ordinului minim
    [n_cheby1, Wn] = cheb1ord(omega_p / pi, omega_s / pi, -20 * log10(1 - Delta_p), -20 * log10(i));
    if n_cheby1~=inf
        ord_cheby1 = [ord_cheby1, n_cheby1];
        [b, a] = cheby1(n_cheby1, -20 * log10(1 - Delta_p), Wn);

        [scor_cost, scor_calitate, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, i);
        scor_cost_cheby1 = [scor_cost_cheby1, scor_cost];
        scor_calitate_cheby1 = [scor_calitate_cheby1, scor_calitate];
        scor_compus_cheby1 = [scor_compus_cheby1, scor_compus];

        delta_s_cheby1 = [delta_s_cheby1, i];

        if scor_compus > scor_max_cheby1(1)
            scor_max_cheby1 = [max(scor_max_cheby1(1), scor_compus), i];
        end
    end

    % folosesc funtia data de MATLAB pentru gasirea ordinului minim
    [n_cheby2, Wn] = cheb2ord(omega_p / pi, omega_s / pi, -20 * log10(1 - Delta_p), -20 * log10(i));
    if n_cheby2~=inf
        ord_cheby2 = [ord_cheby2, n_cheby2];
        [b, a] = cheby2(n_cheby2, -20 * log10(i), Wn);

        [scor_cost, scor_calitate, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, i);
        scor_cost_cheby2 = [scor_cost_cheby2, scor_cost];
        scor_calitate_cheby2 = [scor_calitate_cheby2, scor_calitate];
        scor_compus_cheby2 = [scor_compus_cheby2, scor_compus];

        delta_s_cheby2 = [delta_s_cheby2, i];

        if scor_compus > scor_max_cheby2(1)
            scor_max_cheby2 = [max(scor_max_cheby2(1), scor_compus), i];
        end
    end
end

figure;
scatter(scor_cost_but, scor_calitate_but, 15, '.'); hold on;
scatter(scor_cost_elip, scor_calitate_elip, 15, '.');
scatter(scor_cost_cheby1, scor_calitate_cheby1, 15, '.');
scatter(scor_cost_cheby2, scor_calitate_cheby2, 15, '.');
xlim([0 0.5]); ylim([0 0.5]);
grid on;
legend("But", "Ellip", "Cheby1", "Cheby2");
ylabel("Scor Calitate");
xlabel("Scor Cost");

% din graficele scorului, si din maximele salvate, alegem Delta_s = 0.00749
% aceste este punctul in care unul dintre filtre atinge scorul maxim din toate incercarile
% din for-ul de mai sus reiese Delta_s pentru scorul maxim este approx 0.0075, dar in valoarea fixa, Butterworth iese din tolerante
% asa ca am ales o valoare cu foarte putin mai mica

delta_s_gasit = 0.00749;

b_rezultat = cell(1,4);
a_rezultat = cell(1,4);
titluri    = strings(1,4);

[b, a] = But_FTI(omega_p / pi, omega_s / pi, Delta_p, delta_s_gasit, Ts);
[~, ~, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, delta_s_gasit);
[H1, ~] = freqz(b, a, 5000);
M1 = size(a, 2) - 1; 
H1_name = "Butterworth";
b_rezultat{1} = b;
a_rezultat{1} = a;
titluri(1) = "Filtru Butterworth, $M=" + (size(a, 2) - 1) + "$, $scor=" + scor_compus + "$";

[n_cheby2, Wn] = cheb2ord(omega_p / pi, omega_s / pi, -20 * log10(1 - Delta_p), -20 * log10(delta_s_gasit));
[b, a] = cheby2(n_cheby2, -20 * log10(delta_s_gasit), Wn);
[~, ~, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, delta_s_gasit);
[H2, ~] = freqz(b, a, 5000);
M2 = n_cheby2;
H2_name = "Cheby2";
b_rezultat{2} = b;
a_rezultat{2} = a;
titluri(2) = "Filtru Chebyshev2, $M=" + n_cheby2 + "$ $scor=" + scor_compus + "$";

[n_cheby1, Wn] = cheb1ord(omega_p / pi, omega_s / pi, -20 * log10(1 - Delta_p), -20 * log10(delta_s_gasit));
[b, a] = cheby1(n_cheby1, -20 * log10(1 - Delta_p), Wn);
[~, ~, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, delta_s_gasit);
[H3, ~] = freqz(b, a, 5000);
M3 = n_cheby1; 
H3_name = "Cheby1";
b_rezultat{3} = b;
a_rezultat{3} = a;
titluri(3) = "Filtru Chebushev1, $M=" + n_cheby1 + "$ $scor=" + scor_compus + "$";

[n_elip, Wn] = ellipord(omega_p / pi, omega_s / pi, -20 * log10(1 - Delta_p), -20 * log10(delta_s_gasit));
[b, a] = ellip(n_elip, -20 * log10(1 - Delta_p), -20 * log10(delta_s_gasit), Wn);
[scor_cost, scor_calitate, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, delta_s_gasit);
[H4, ~] = freqz(b, a, 5000);
M4 = n_elip;
H4_name = "Cauer"; 
b_rezultat{4} = b;
a_rezultat{4} = a;
titluri(4) = "Filtru Cauer, $M=" + n_elip + "$ $scor=" + scor_compus + "$";

save("faza3.mat", "H1", "H2", "H3", "H4", "M1", "M2", "M3", "M4", "H1_name", "H2_name", "H3_name", "H4_name");

%% Plot

function grafic_faza_3(b, a, omega_p, omega_s, Delta_p, Delta_s, rez, titlu)

% functie locala pentru a a desena matricea de functii de transfer ale celor 4 filtre

%% flag-uri globale

global EXTRA;
global DEBUG;

%% Corpul functiei
figure;

for i = 1:4
    subplot(2, 4, i);
    [H, w] = freqz(b{i}, a{i}, rez);
    faza = unwrap(angle(H));
    valori_corecte = abs(H) > eps;

    %% Plot caracteristica frecventa
    subplot(2, 4, i);

    H_abs=abs(H);
    plot(w, db(H_abs));

    linii_verticale = [omega_p, omega_s];
    linii_orizontale = [1 + Delta_p, 1 - Delta_p, Delta_s];

    xline(linii_verticale);
    yline(db(linii_orizontale));

    ylabel("$|G(e^{\omega j})|,\ in\ dB$", 'interpreter', 'latex');
    xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
    xlim([0, pi]);
    ylim([db(eps), 10]);
    title(titlu(i), 'Interpreter','latex');

    %% Plot caracteristica in faza

    subplot(2, 4, i + 4);
    plot(w(valori_corecte), faza(valori_corecte));

    xline(linii_verticale);

    ylabel("$Faza\ (\ radian\ )$", 'interpreter', 'latex');
    xlabel("$\omega\ (\ radian\ )$", 'interpreter', 'latex');
    xlim([0, pi]);

end
end

grafic_faza_3(b_rezultat, a_rezultat, omega_p, omega_s, Delta_p, delta_s_gasit, 5000, titluri);