% Balas Vlad-George, 333AB

omega_p = 0.9595
omega_s = 1.2633
Delta_p = 0.0623
Delta_s = Delta_p % in conformitate cu indicatiie date
Ts = 1.7027

rez = 5000

%% Functii utilitare
tustin = @(omega, t_s)((2 / t_s) * tan(omega / 2));
inverse_tustin = @(ohm, t_s)( 2 * atan((t_s / 2) * ohm));

pseudo_tustin = @(omega, t_s)((1 / t_s) * tan(omega / 2));
inverse_pseudo_tustin = @(ohm, t_s)(2 * atan(t_s * ohm));

%% Faza 4

[B_faza_4, A_faza_4] = fnc(omega_p/ pi, omega_s / pi, Delta_p, Delta_s, Ts, "faza_4");
M = size(A_faza_4, 2) - 1;
M_p = 1 - Delta_p;

ohm_p = tustin(omega_p, Ts);
ohm_c = ohm_p / nthroot((Delta_p + 2) * Delta_p, 2 * M);
omega_c = inverse_tustin(ohm_c, Ts);

% matrice linii verticale
xlines = [omega_p / pi, omega_s / pi , omega_c / pi;
    omega_p / pi, omega_s / pi, omega_c / pi];

% matrice linii orizontale
ylines = [1 + Delta_p, 1 - Delta_p, Delta_s;];

grafic_abs(B_faza_4, A_faza_4, xlines, ylines, "Caracteristica filtrului pentru Faza 4");