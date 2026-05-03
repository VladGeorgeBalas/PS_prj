% Balas Vlad-George, 333AB

% Salvate dupa rularea functiei
omega_p = 0.9595
omega_s = 1.2633
Delta_p = 0.0623
Delta_s = 2 * Delta_p % in conformitate cu indicatiie date
Ts = 1.7027

tustin = @(omega, t_s)((2 / t_s) * tan(omega / 2));
inverse_tustin = @(ohm, t_s)(2 * atan((t_s / 2) * ohm));
M_p = 1 - Delta_p;

[B, A] = fnc(omega_p / pi, omega_s / pi, Delta_p, Delta_s, Ts, "tustin");
M = size(A, 2) - 1;
M_dat = M;
omega_c = inverse_tustin(...
    tustin(omega_p, Ts) / nthroot((1 - M_p * M_p)/(M_p * M_p), 2 * M),...
    Ts);

%% Subpunctul a
disp("elliptic");
[M, b, a] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M_dat, "ellip");
grafic_abs(b, a, [omega_p/pi, omega_c/pi, omega_s/pi; omega_p/pi, omega_c/pi, omega_s/pi], [1 + Delta_p, 1 - Delta_p, Delta_s], "Filtul eliptic de ordin M = " + sprintf('%d', M));

%% Subpuctul b

disp("Butterworth")
grafic_abs(B, A, [omega_p/pi, omega_c/pi, omega_s/pi; omega_p/pi, omega_c/pi, omega_s/pi], [1 + Delta_p, 1 - Delta_p, Delta_s], "Filtul Butterworth de ordin M = " + sprintf('%d', M_dat));
[scor_cost_mediu, scor_tol_mediu, scor_compus] = criteriu(B, A, omega_p, omega_s, Delta_p, Delta_s)

disp("fir1");
[M, b, a] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M_dat, "fir1");
grafic_abs(b, a, [omega_p/pi, omega_c/pi, omega_s/pi; omega_p/pi, omega_c/pi, omega_s/pi], [1 + Delta_p, 1 - Delta_p, Delta_s], "Filtul FIR folosind metoda ferestrei de ordin M = " + sprintf('%d', M));
[scor_cost_mediu, scor_tol_mediu, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, Delta_s)

disp("firls");
[M, b, a] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M_dat, "firls");
grafic_abs(b, a, [omega_p/pi, omega_c/pi, omega_s/pi; omega_p/pi, omega_c/pi, omega_s/pi], [1 + Delta_p, 1 - Delta_p, Delta_s], "Filtul FIR folosint metoda CMMP de ordin M = " + sprintf('%d', M));
[scor_cost_mediu, scor_tol_mediu, scor_compus] = criteriu(b, a, omega_p, omega_s, Delta_p, Delta_s)

%% Subpunctul c

disp("cheby1");
[M, b, a] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M_dat, "cheby1");
grafic_abs(b, a, [omega_p/pi, omega_c/pi, omega_s/pi; omega_p/pi, omega_c/pi, omega_s/pi], [1 + Delta_p, 1 - Delta_p, Delta_s], "Filtul Chebisev de ordin M = " + sprintf('%d', M));

disp("cheby2");
[M, b, a] = bin_opt(Delta_p, Delta_s, omega_p, omega_s, Ts, M_dat, "cheby2");
grafic_abs(b, a, [omega_p/pi, omega_c/pi, omega_s/pi; omega_p/pi, omega_c/pi, omega_s/pi], [1 + Delta_p, 1 - Delta_p, Delta_s], "Filtul Chebisev de ordin M = " + sprintf('%d', M));
