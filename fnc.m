% Balas Vlad-George, 333AB
%--------------------------
% FILE: fnc.m 
% 
% FUNCTION: fnc
%
% CALL: [B, A]	= fnc(omega_p, omega_s, delta_p, delta_s, t_s, metoda)
% 
% Functia returneaza numitorul si numaratorul functiei de transfer al unui filtru IIR
% Butterworth obtinut cu specificatiile date si metoda de transformare ceruta.
%
% Argumente
%	omega_p	- Frecventa banda de trecere, in [0, 1]
%	omega_s	- Frecventa banda de stopare, in [0, 1]
%	delta_p	- Toleranta banda de trecere, [0, 1]
%	delta_s	- Toleranta banda de stopare, [0, 1]
%	t_s		- Perioada de esantionare
%	metoda	- metoda prin care se obtine filtrul: 'tustin', 'pseudo_tustin' sau 'faza_4'
%
% Iesire
%	B		- Numaratorul functiei de transfer
%	A		- Numitorul functiei de transfer
%
% Foloseste
%	WAR_ERR
%
% Autor: Balas Vlad-George
% Creat: Decembrie, 2025

function [B,A] = fnc(omega_p,omega_s,delta_p,delta_s,t_s, metoda)
FN = '<BUT_FTI>: ' ;
E1 = [FN 'Missing, empty or inconsistent input data => empty outputs. Exit.'] ;

% Flag-uri globale pentru a putea seta modul de rulare al codului
global EXTRA					% flag-ul extra il folosesc sa printez informatii extra
global DEBUG					% flag-ul debug printeaza mai multe informatii pentru cand nu functioneaza

%% Verificarea corectitudinii argumentelor
B = [];
A = [];

if (nargin < 3)
	war_err(E1);
	return;
end

omega_p = abs(omega_p(1)) ;
if (omega_p < eps) || (omega_p >= (1-eps))
	war_err(E1);
	return;
end

omega_s = abs(omega_s(1)) ;
if (omega_s < eps) || (omega_s >= (1-eps))
	war_err(E1);
	return;
end

delta_p = abs(delta_p(1)) ;
if (delta_p < eps) || (delta_p >= (1-eps))
	war_err(E1);
	return;
end

if (nargin < 4)
	delta_s = delta_p ;
end

delta_s = abs(delta_s(1)) ;
if (delta_s < eps) || (delta_s >= (1-eps))
	war_err(E1);
	return;
end

if (omega_p > omega_s)
	FN = omega_p ;
	omega_p = omega_s ;
	omega_s = FN ;
end

if (nargin < 5)
	t_s = 2 ;
end

t_s = abs(t_s(1)) ;
if (t_s < eps)
	t_s = 2 ;
end

% Verific daca e o metoda implementata si daca nu, rezolva pt tustin
if ~(strcmp(metoda, "tustin") == 1 || strcmp(metoda, "pseudo_tustin") == 1 || strcmp(metoda, "faza_4") == 1)
    metoda = "tustin";
end

%% Corpul functiei

if (strcmp(metoda, "tustin") == 1)
	
	% Acesta este corpul functiei BUT_FTI. Am notat variabilele cum erau
	% indicate in comentarii. In rest, nu am modificat nimic.

	% Am vrut sa o am aici ca sa fie clar unde si ce modificari am facut
	% pentru rezolvarile de la 1e si 4


	if DEBUG
		disp("tustin")
	end

	% Functii utilitare
	tustin = @(omega, t_s)((2 / t_s) * tan(omega * pi / 2));
	inverse_tustin = @(ohm, t_s)((2 / pi) * atan((t_s / 2) * ohm));

	ohm_p = tustin(omega_p, t_s); 										% Compute Omega_p.
	ohm_s = tustin(omega_s, t_s); 										% Compute Omega_s.
	m_p = 1-delta_p;

	m_p_pow2 = m_p * m_p;
	delta_s_pow2 = delta_s * delta_s;

	M = ceil(log((m_p_pow2 * (1 - delta_s_pow2))....
		/(delta_s_pow2 * (1 - m_p_pow2)))...
		/ (2 * log(ohm_s/ohm_p)));

	ohm_c = ohm_p / nthroot((1 - m_p_pow2)/m_p_pow2, 2 * M);

	s_m = ohm_c * exp(1j*(M+(1:2:(2*M)))*pi/(2*M)); 					% The stable poles of filter.
	s_m = s_m * t_s ;

	FN = s_m;															% pentru compatibilitate de notatie
	E1 = 2-FN ;
	B = real(prod(-FN./E1)*poly(-ones(1,M))) ; 							% Numerator of transfer function.
	A = real(poly((2+FN)./E1)) ; 										% Denominator of transfer function.
elseif (strcmp(metoda, "pseudo_tustin") == 1)
	if DEBUG
		disp("pseudo_tustin")
	end

	% Functii utilitare
	pseudo_tustin = @(omega, t_s)((1 / t_s) * tan(omega * pi / 2));
	inverse_pseudo_tustin = @(ohm, t_s)((2 / pi) * atan(t_s * ohm));

	ohm_p = pseudo_tustin(omega_p, t_s); 	% Compute Omega_p.
	ohm_s = pseudo_tustin(omega_s, t_s); 	% Compute Omega_s.
	m_p = 1 - delta_p;

	m_p_pow2 = m_p * m_p;
	delta_s_pow2 = delta_s * delta_s;

	M = ceil(log((m_p_pow2 * (1 - delta_s_pow2))/(delta_s_pow2 * (1 - m_p_pow2)))...
		/ (2 * log(ohm_s/ohm_p)));
	ohm_c = ohm_p / nthroot((1 - m_p_pow2)/m_p_pow2, 2 * M);

	s_m = ohm_c * exp(1j*(M+(1:2:(2*M)))*pi/(2*M));
	s_m = s_m * t_s ;

	B = real(prod(-s_m./(1 - s_m))*poly(-ones(1,M))) ; 	% Numerator of transfer function.
	A = real(poly((1 + s_m)./(1 - s_m))) ; 			% Denominator of transfer function.
elseif (strcmp(metoda, "faza_4") == 1)
	if DEBUG
		disp("faza_4")
	end

	% Functii utilitare
	tustin = @(omega, t_s)((2 / t_s) * tan(omega * pi / 2));
	inverse_tustin = @(ohm, t_s)((2 / pi) * atan((t_s / 2) * ohm));

	ohm_p = tustin(omega_p, t_s); 										% Compute Omega_p.
	ohm_s = tustin(omega_s, t_s); 										% Compute Omega_s.
	
	M = ceil(log((1 + delta_p - delta_s) * (1 + delta_p + delta_s)...	% Ordinul functiei de transfer
	 	/ (delta_p * delta_s * delta_s))...
		/ (2 * log(ohm_s/ohm_p)));

	ohm_c = ohm_p / nthroot((delta_p + 2) * delta_p, 2 * M);			% Omega_c conform formulei gasite

	s_m = ohm_c * exp(1j*(M+(1:2:(2*M)))*pi/(2*M));						% Gasirea polilor
	s_m = s_m * t_s ;

	B = (1 + delta_p) * real(prod(-s_m./(2 - s_m))*poly(-ones(1,M))) ; 	% Numerator of transfer function.
	A = real(poly((2 + s_m)./(2 - s_m))) ; 								% Denominator of transfer function.

end
end