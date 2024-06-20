clear;
clc;
close all;
%% Simulation parameters
true_DoAs = [0 -20 50];  % True DoA angles (in degrees)
%For better accuracy in FFT method, increase no. of antennas M. MUSIC gives
%accurate results for any number of antennas
M = 8;                        
K = length(true_DoAs);        % Number of sources 
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
    x = A; % * S + n; % Received signal matrix

   theta = -90:0.1:90; % ULA angle sweep
   Rxx = (x * x') / size(x, 1); % Covariance Matrix 

 
   [theta_fft, fft_peaks, fft_estimates, Yq] = InterpFFT(Rxx, M, K);
   [doa_peaks, doa_estimates, error_function, Pmusic] = MUSIC(Rxx, K, d, M, wavelength, theta);

    %% Plots

    plot(theta_fft, abs(Yq));
    xlabel('Angle (degrees)');
    ylabel('FFT Spectrum');
    title('Interpolated FFT');

    figure;
    fplot(error_function); 
    title("Error function b/w actual DoAs and estimated MUSIC DoAs");

    figure;
    plot(theta, 10 * log10(abs(Pmusic)/max(abs(Pmusic))), 'LineWidth', 2);
    xlabel('Angle (degrees)');
    ylabel('Normalized Power (dB)');
    title('MUSIC Spectrum');
    grid on;

    % Plot the true and estimated DoAs
    figure;
    plot(1:K, sort(true_DoAs), 'ro', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'True DoAs');
    hold on;
    plot(1:K, sort(doa_estimates), 'bx', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Estimated DoAs (MUSIC)');
    hold on;
    plot(1:K, sort(fft_estimates), 'diamond', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'Estimated DoAs (InterpFFT)');
    xlabel('Source index');
    ylabel('Angle (degrees)');
    title('True and Estimated DoAs');
    legend;
    grid on;
    hold off;
% Error between actual and estimated DoAs
    error_MUSIC = sort(true_DoAs) - sort(doa_estimates);
    display("Error between true and MUSIC DoAs:", num2str(error_MUSIC));
    error_FFT = sort(true_DoAs) - sort(fft_estimates);
    display("Error between true and MUSIC DoAs:", num2str(error_FFT));

  %% Interpolated FFT
  function [theta_fft, fft_peaks, fft_estimates, Yq] = InterpFFT(Rxx, M, K)
   Y = fftshift(fft(Rxx)); 
   Yq = zeros(1, length(0:0.05:length(Y(:,1))-1));
   for i = 1 : M
   Yq = Yq + interp1(0:length(Y(:,i))-1, abs(Y(:, i)), 0:0.05:length(Y(:,i))-1, 'cubic');
   end
   Yq = fliplr(Yq);  
   theta_fft=-90:180/length(Yq):90-1/length(Yq);
   figure;
   [fft_peaks, fft_estimates] = findpeaks(Yq, theta_fft);
    [~, indexfft] = sort(fft_peaks, 'descend');
    fft_estimates = fft_estimates(:, indexfft);
    if(length(fft_estimates) > K)
     fft_estimates = fft_estimates(1:K);
    end
    
     % If less than K estimates were found, fill the remaining ones with NaN
    if length(fft_estimates) < K
        fft_estimates = [fft_estimates, nan(1, K - length(fft_estimates))];
    end 
  end

%% MUSIC
function [doa_peaks, doa_estimates, error_function, Pmusic] = MUSIC(Rxx, K, d, M, wavelength, theta)
   [V, D] = eig(Rxx); % Eigen decomposition

   % Arrange eigen values in descending order - The eigen vectors corresponding to the first 'K' eigenvalues will
   % represent signal (+ low SNR noise) subspace while the remaining
   % (M-K) eigenvectors will comprise the noise subspace
   [~, idx] = sort(diag(D), 'descend'); % Arrange eigenvectors in descending order
   V = V(:, idx); % Arrange eigenvectors corresponding to the order of the eigenvalues above
   Vn = V(:, K+1:end); % Pure noise subspace
   

    % Calculate the MUSIC spectrum by implementing the pseudospectrum
    % equation
    Pmusic = arrayfun(@(angle) 1 / abs((a(angle, d, M, wavelength)' * Vn) * (Vn' * a(angle, d, M, wavelength))), theta);
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
    [doa_peaks, doa_estimates] = findpeaks(Pmusic, theta);
    [~, indexpeak] = sort(doa_peaks, 'descend');
    doa_estimates = doa_estimates(:, indexpeak);
    if(length(doa_estimates) > K)
     doa_estimates = doa_estimates(1:K);
    end
    % If less than K estimates were found, fill the remaining ones with NaN
    if length(doa_estimates) < K
        doa_estimates = [doa_estimates, nan(1, K - length(doa_estimates))];
    end
end

%% Function to compute steering vector
function a_vec = a(angle, d, M, wavelength)
    a_vec = exp(-1i*2*pi*d*(0:M-1)'*sind(angle)/wavelength);
end

