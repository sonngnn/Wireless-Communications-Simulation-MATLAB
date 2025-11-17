clc;
close all;
clear all;

% Paramètres
Ns = 5000;
Fe = 4*10^6;
Fse = 4;
span = 8;
alpha = 0.35;
Nfft = 1024;

% Constellation d'un QPSK en mapping de Gray
constellation_4 = mod_psk(4, true);
figure;
plot_cstl(constellation_4);
title('Constellation QPSK');

% ------------------------ EMETTEUR ------------------------
mapping_x = randi([0 3],1,Ns);
x = constellation_4(mapping_x+1);

% Suréchantillonage et mise en forme
x_upsampled = upsample(x, Fse);
g = rcosdesign(alpha,span,Fse);
sl = conv(x_upsampled, g);

% Définition des paramètres des différents canaux à tester
canaux = [
    0, 1, 0;    % (α = 0, d = 1, ϕ = 0)
    1, 1, 0;    % (α = 1, d = 1, ϕ = 0)
    1, 4, 0;    % (α = 1, d = 4, ϕ = 0)
    1, 1, pi;   % (α = 1, d = 1, ϕ = π)
    0.5, 1, 0;  % (α = 0.5, d = 1, ϕ = 0)
    0.5, 4, 0   % (α = 0.5, d = 4, ϕ = 0)
];

colors = ['r', 'g', 'b', 'y', 'c', 'm'];

% ------------------------ CANAL ------------------------
figure;
hold on;
title('Constellations après filtrage adapté');

for i = 1:size(canaux,1)
    alpha_c = canaux(i,1);
    d_c = canaux(i,2);
    phi_c = canaux(i,3);
    
    % Canal à 2 trajets
    h_canal = zeros(1, d_c + 1);
    h_canal(1) = 1;
    h_canal(d_c + 1) = alpha_c * exp(1i * phi_c);
    
    % Normalisation
    C = sum(abs(h_canal).^2);
    h_canal = h_canal / C;
    
    % Passage du signal à travers le canal
    y = conv(sl, h_canal);   

    mu = 0;
    sigma=0;
    zn = randn(1,length(y));
    zn = mu + sigma*zn;

    y = zn + y;

    % ------------------------ RECEPTEUR ------------------------
    
    % Filtrage adapté
    g_adapte = fliplr(g);
    rl = conv(y, g_adapte);
    
    % Échantillonnage
    indice_echantillonnage = length(g);
    rl_echant = rl(indice_echantillonnage:Fse:indice_echantillonnage+(Ns-1)*Fse);
    
    % Tracé de constellation après réception
    subplot(2,3,i);
    plot(rl_echant, '.', 'Color', colors(i), 'MarkerSize', 10);
    grid on;
    title(['α=', num2str(alpha_c), ', d=', num2str(d_c), ', ϕ=', num2str(phi_c)]);
    xlabel('Re');
    ylabel('Im');
    xlim([-2 2]);
    ylim([-2 2]);
end
hold off;


% ------------------------ TRAÇAGE DSP ------------------------
f = linspace(-Fe/2, Fe/2, Nfft); % Axe des fréquences

figure,
hold on;
title('DSP du signal transmis et réponses fréquentielles des canaux');
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
grid on;

% Calcul et affichage de la DSP du signal transmis avec pwelch
[pxx, f_pxx] = pwelch(sl, hamming(Nfft), [], Nfft, Fe, 'centered');
plot(f_pxx, pxx*Fe, 'k', 'LineWidth', 1, 'DisplayName', 'DSP du signal transmis');

for i = 1:size(canaux,1)
    alpha_c = canaux(i,1);
    d_c = canaux(i,2);
    phi_c = canaux(i,3);
    
    % Canal à 2 trajets
    h_canal = zeros(1, d_c + 1);
    h_canal(1) = 1;
    h_canal(d_c + 1) = alpha_c * exp(1i * phi_c);
    
    % Normalisation
    C = sum(abs(h_canal).^2);
    h_canal = h_canal / C;
    
    % Réponse fréquentielle du canal (FFT)
    H_f = fftshift(abs(fft(h_canal, Nfft))); % Magnitude de la réponse fréquentielle

    % Tracé des réponses fréquentielles des canaux
    plot(f, H_f, 'Color', colors(i), 'LineWidth', 1, 'DisplayName', ['α=', num2str(alpha_c), ', d=', num2str(d_c), ', ϕ=', num2str(phi_c)]);
end

legend show;
hold off;

% ------------------------ EGALISATION ------------------------
v = conv(g_adapte, conv(h_canal, g));
indice_echantillonnage = length(g);
v_echant = v(indice_echantillonnage:Fse:end);
figure;
plot(v_echant);

L = length(v_echant); 
N = 2*L;
V = zeros(N, N + L - 1);

for i = 1:N
    V(i, i:i+L-1) = flip(v_echant);
end
