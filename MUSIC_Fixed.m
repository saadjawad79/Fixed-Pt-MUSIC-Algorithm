
clear;
clc;
close all;
%% Simulation parameters
true_DoAs = [0, -20 50];  % True DoA angles (in degrees)
M = 8;                        % Number of antennas (fixed at 8 for LUTs)
K = length(true_DoAs);        % Number of sources (fixed at 3)
c = 3e9;                      % Speed of light
f = 28e9;                     % Frequency
wavelength = c/f;             % Wavelength(lambda)
d = wavelength / 2;           % Distance between antennas (half-wavelength)
SNR = 20;                     % Signal-to-noise ratio (dB)

% Check if M >= K
if M < K
    error('Number of antennas must be greater than or equal to number of sources');
end

%% Signal vector data generation 
    N = 1000; % Number of samples
    A = exp(-1i*2*pi*d*(0:M-1)'*sind(true_DoAs)/wavelength); % Steering matrix
    S = sqrt(0.5) * (randn(K, N) + 1i * randn(K, N)); % Complex Gaussian signals (can be sinusoids)
    noise_variance = 10^(-SNR/10) * trace(A * (S * S') * A') / (M * N);
    n = sqrt(noise_variance / 2) * (randn(M, N) + 1i * randn(M, N)); % Complex Gaussian noise
    x = A * S + n; % Received signal matrix

M = fi(8,0, 8, 0);                        % Number of antennas (increase for better accuracy in InterpFFT, MUSIC is super resolution so no need)
K = fi(length(true_DoAs), 0, 8, 0);        % Number of sources 
c = 3e9;                      % Speed of light
f = 28e9;                     % Frequency
wavelength = fi(c/f, 0, 32, 18);             % Wavelength(lambda)
d = fi(wavelength/2, 0, 32, 18);           % Distance between antennas (half-wavelength)
SNR = 20;                     % Signal-to-noise ratio (dB)

% Check if M >= K
if M < K
    error('Number of antennas must be greater than or equal to number of sources');
end

%% Signal vector data generation 
    N = 1000; % Number of samples
    A = fi(A, 1, 32, 14); % Steering matrix
    S = fi(S, 1, 32, 14); % Complex Gaussian signals (can be sinusoids)
    noise_variance = fi(noise_variance, 1, 32, 14);
    n = fi(n, 1, 32, 14); % Complex Gaussian noise
    x = fi(x, 1, 32, 14); % Received signal matrix

   theta = fi(-90:0.1:90, 1, 32, 14); % ULA angle sweep
   Rx = (x * x') / size(x, 1); % Covariance Matrix 
   Rxx = fi(Rx, 1, 32, 14);
  %% MUSIC
   [V, D] = eig(double(Rxx)); % Eigen decomposition
   V = fi(V, 1, 32, 14);
   D = fi(D, 1, 32, 14);
 
   [~, idx] = sort(diag(D), 'descend'); % Arrange eigenvectors in descending order
   V = V(:, idx); % Arrange eigenvectors corresponding to the order of the eigenvalues above
   Vn = V(:, K+1:end); % Pure noise subspace
   

    % Calculate the MUSIC spectrum by implementing the pseudospectrum
    % equation
    Pmusic = fi(zeros(1,length(theta)), 1, 64, 50);
    Pmusic(:) = arrayfun(@(angle) 1 / abs((a(angle, d, M, wavelength)' * Vn) * (Vn' * a(angle, d, M, wavelength))), theta);
   
    % Find the initial estimates of DoAs
    [~, doa_indices] = sort(Pmusic, 'descend');

    % If more than K peaks are found, select the K largest ones
    if length(doa_indices) > K
        doa_indices = doa_indices(1:K);
    end

    initial_DoAs = theta(doa_indices);

    for i = 1:length(initial_DoAs)
        error_function = @(x) abs(1 / abs((a(x, d, M, wavelength)' * Vn) * (Vn' * a(x, d, M, wavelength))) - Pmusic(doa_indices(i)));
    end
 
    
    % Find DoA estimates for MUSIC
    [doa_peaks, doa_estimates] = findpeaks(double(Pmusic), double(theta));
    [~, indexpeak] = sort(doa_peaks, 'descend');
    doa_estimates = doa_estimates(:, indexpeak);
    inv_max = 1/max(Pmusic);
    if(length(doa_estimates) > K)
     doa_estimates = doa_estimates(1:K);
    end
    % If less than K estimates were found, fill the remaining ones with NaN
    if length(doa_estimates) < K
        doa_estimates = [doa_estimates, nan(1, K - length(doa_estimates))];
    end
    % Plot (only for simulation purposes)
    figure;
    plot(theta, Pmusic*inv_max);
    figure;
    dBplot = zeros(1, length(theta));
    for i = 1: length(theta)
    dBplot(i) = fi(replacement_log10((Pmusic(i))*inv_max), 1, 32, 16);
    end
    plot(theta, 10 * dBplot, 'LineWidth', 2);
    xlabel('Angle (degrees)');
    ylabel('Normalized Power (dB)');
    title('MUSIC Spectrum');
    grid on;

    
    figure;
    plot(1:K, sort(true_DoAs), 'ro', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'True DoAs');
    hold on;
    plot(1:K, sort(doa_estimates), 'bx', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Estimated DoAs (MUSIC)');
    hold on;
    xlabel('Source index');
    ylabel('Angle (degrees)');
    title('True and Estimated DoAs');
    legend;
    grid on;
    hold off;
    % Error between actual and estimated DoAs
    error_doas = sort(true_DoAs) - sort(doa_estimates);
    display(num2str(error_doas));

%% Function to compute steering vector
function a_vec = a(angle, d, M, wavelength)
    a_vec = cordiccexp(mod(-(0:M-1)'*fi(pi, 1, 32, 14)*replacement_sind(angle), fi(2*pi, 1, 32, 14)));
    % a_vec = cos(-pi*(0:double(M)-1)'*double(replacement_sind(angle))) + 1i*sin(-pi*(0:double(M)-1)'*double(replacement_sind(angle))); 
    % a_vec = exp(-1i*2*pi*double(d)*(0:double(M)-1)'*double(replacement_sind(angle))/double(wavelength));
end

