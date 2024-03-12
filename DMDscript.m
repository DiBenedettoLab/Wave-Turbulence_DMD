% This MATLAB script performs wave-turbulence decomposition on raw flow
% measurement data. It includes spectral analysis, Proper Orthogonal
% Decomposition (POD), Dynamic Mode Decomposition (DMD), and energy
% distribution checks. The script visualizes spatial modes, DMD spectra,
% and reconstructs time series to separate wave and turbulence components.

% Input:
%   - Raw flow measurement data loaded from 'data.mat'
%   - Sampling frequency, fs
%   - Wave frequency range
%   - For spectra: number of fft points, nfft and window size, win

% Output:
%   - Various plots showcasing PSDs, spatial modes, DMD spectrum, and energy
%     distribution.

% Tasks:
% 1. Load raw data and set up parameters.
% 2. Visualize power spectra and adjust wave frequency range.
% 3. Construct Hankel matrix and apply POD to analyze fluctuation data 
%    to determine rank truncation r.
% 4. Perform DMD to obtain modes and eigenvalues and find modes in 
%    wave frequency range.
% 5. Check energy distribution through eigenvalues and frequencies.

% Dependencies
% This script needs MATLAB signal processing toolbox 23.2 or later
%   - HankelMatrix: Function for constructing Hankel matrix.
%   - DMD_processing: Function for processing raw input signal and performing DMD.
%       - DMDselective: Function for DMD analysis.
%       - HankelMatrix: Hankel matrix construction.
%   - DMD_spectrum: Function for computing the DMD spectrum.
%   - find_in_range: Function for finding indexes of values within wave range.
%   - DMD_reconstruct: Function for reconstructing signals using DMD.
%       - textprogressbar: Utility function for displaying progress bars.
%   - eigvals_lambda_omega: Function for extracting components and angles from eigenvalues.


%% Load data
close
clear
clc
%==========================================================================
nfft = 2048; win = 1000;
load('SynthData.mat', 'u_tot'); % 1xk time series of flow measurements (ideally in m/s)
DataRaw = u_tot;
fs = 10; % sampling frequency (Hz)
waveRange = [1.17e-1 6e-1]; % [low-freq_end hig-freq_end]
%==========================================================================
dt = 1/fs; 
[numSignals, numObs] = size(DataRaw);
t = (0:numObs-1)/fs;


%% plot spectra of fluctuation and adjust waveRange accordingly

DataFluc = DataRaw - mean(DataRaw);
% spectra of fluctuations
[Sxx, fm] = cpsd(DataFluc', DataFluc',hann(win),win/2,nfft,fs);

figure(4);
loglog(fm, Sxx, 'Color', 'k', 'LineWidth',1)
hold off
xlim([fm(1) fm(end)])
grid on
xlabel('Frequency (Hz)')
ylabel('PSD (m$^2$s$^{-2}$/Hz)')

% section of time series
figure(5);
ts = 500; ss = 1500;
plot(t(1,ts:ts+ss),DataFluc(1,ts:ts+ss), 'Color', 'k', 'LineWidth',1)
xlabel('Time (s)', 'Interpreter','latex')
ylabel('Velocity (m/s)', 'Interpreter','latex')
xlim([t(ts) t(ts+ss)])


%% Decomposition
clc
% ========================================================================
rawInputSignal = DataFluc; % Should not have mean
r =[1 131]; % range of modes to use in truncation
matParam.Lag = 1; % number of obs to lag
matParam.NumObsLag = 1500; % number of columns in time-delayed matrix
% ========================================================================

% POD spectrum

[X1, ~] = HankelMatrix(rawInputSignal,matParam.NumObsLag ,numObs - matParam.NumObsLag, 'column');
% Decompose time-delayed matrix using SVD
[U, S, V] = svd(X1, 'econ');
V_win = round(size(V, 1)/3);
Vspec = nan(size(V, 1), nfft/2 + 1);
for ii = 1:size(V, 1)
    % get power spectrum of each row of V (or each column of V')
    Vspec(ii,:) = cpsd(V(:,ii),V(:,ii),hann(V_win),round(V_win/2),nfft,fs);
end
figure(12);
[X,Y] = meshgrid(fm,1:size(V, 1));
surf(X,Y,Vspec, 'EdgeColor', 'none')
% colormap('sky')
colormap(flipud(gray))
set(gca,'xscale','log')
set(gca,'yscale','log')
view(2); % Set the viewing angle to show only the x and y axes
box on
grid off
xline(waveRange, ':')
yline(r,'--')
xlabel('Frequency (Hz)', 'Interpreter','latex')
ylabel('Mode $j$', 'Interpreter','latex')
xlim([fm(1) fm(end)])

%% DMD
[X1, X2, Out] = DMD_processing(rawInputSignal, r, matParam, dt);

% DMD spectrum
[P2, f2, P1, f1] = DMD_spectrum(Out.b, Out.lambda, dt);
% find the indices of values within 'waveRange' in the input vector 'f1' 
% for 1-sided spectra and 'f2' for two-sided spectra.
[in1,out1, ~, ~] = find_in_range(f1,waveRange);
[in2,out2, inPos2, outPos2] = find_in_range(f2,waveRange);

figure(101)
loglog(f1(in1), P1(in1), 'xr'); hold on
loglog(f1(out1), P1(out1), 'xk');
xline(waveRange, '--')
box on
ylabel('DMD spectrum', 'Interpreter','latex'); xlabel('Frequency (Hz)', 'Interpreter','latex')


[u_time_dynamics_wave, u_Xdmd_wave] = DMD_reconstruct(X1,Out.Phi(:,in2),Out.omega(in2), Out.b(in2),dt);

u_wave_dmd = [real(u_Xdmd_wave(1,:)) real(u_Xdmd_wave(2:end,end))'];

u_turb_dmd = rawInputSignal(1:end-1) - u_wave_dmd;

% Plotting wave time series (DMD)
figure(3);
hold on
s1 = 2000; s2 = 2800;
plot(t(1:s2-s1+1), rawInputSignal(1,s1:s2), 'Color', '#444444', 'LineStyle','-', 'LineWidth', 1)
plot(t(1:s2-s1+1), u_wave_dmd(s1:s2), 'Color', '#1982c4', 'LineStyle','-', 'LineWidth', 1)
plot(t(1:s2-s1+1), u_turb_dmd(s1:s2), 'Color', '#ff595e', 'LineStyle','-', 'LineWidth', 1)
hold off
xlim([t(1) t(s2-s1+1)])
ylabel('$u$ (m/s)', 'Interpreter','latex'); xlabel('Time (seconds)', 'Interpreter','latex')
legend('Raw signal', 'DMD wave', 'DMD turbulence', 'Interpreter','latex')
box on


% Main DMD spectra
figure(4);
[SxxXs, fmXs] = cpsd(rawInputSignal, rawInputSignal,hann(win),win/2,nfft,fs);
[SxxXfs, fmXfs] = cpsd(u_wave_dmd, u_wave_dmd,hann(win),win/2,nfft,fs);
[SxxXnfs, fmXnfs] = cpsd(u_turb_dmd, u_turb_dmd,hann(win),win/2,nfft,fs);

loglog(fmXs, SxxXs, 'k', 'LineWidth', 1.5); hold on
loglog(fmXfs, SxxXfs, 'Color', '#1982c4', 'LineStyle','-', 'LineWidth', 1)
loglog(fmXnfs, SxxXnfs, 'Color', '#ff595e', 'LineStyle','-', 'LineWidth', 1)

xline(waveRange, ':')
xlim([fmXs(1) fmXs(end)])
hold off
xlabel('Frequency (Hz)', 'Interpreter','latex')
ylabel('PSD (m$^2$s$^{-2}$/Hz)', 'Interpreter','latex')
legend('Raw','DMD wave', 'DMD turbulence', 'Interpreter','latex')


%% Check eigenvalues
[lambdaCompsIn, lambdaAngleIn, ~, omegaCompsIn] = eigvals_lambda_omega(Out.lambda(in2), Out.omega(in2));
[lambdaCompsOut, lambdaAngleOut, circleArcOut, omegaCompsOut] = eigvals_lambda_omega(Out.lambda(out2), Out.omega(out2));

[~, lambdaAngleInPos, ~, ~] = eigvals_lambda_omega(Out.lambda(inPos2), Out.omega(inPos2));
[~, lambdaAngleOutPos, ~, ~] = eigvals_lambda_omega(Out.lambda(outPos2), Out.omega(outPos2));

figure(10);
subplot(1,2,1)
plot(circleArcOut(1,:),circleArcOut(2,:),'k') % plot unit circle
hold on, grid on
scatter(lambdaCompsIn(1,:),lambdaCompsIn(2,:),'or'); hold on
scatter(lambdaCompsOut(1,:),lambdaCompsOut(2,:),'ok')
xlabel('Real($\lambda$)', 'Interpreter','latex'); ylabel('Imaginary($\lambda$)', 'Interpreter','latex')
axis equal
subplot(1,2,2)
hold on, grid on
scatter(omegaCompsIn(1,:),omegaCompsIn(2,:),'+r')
scatter(omegaCompsOut(1,:),omegaCompsOut(2,:),'+k')
xline(0); yline(0)
box on
xlabel('Real($\omega$)', 'Interpreter','latex'); ylabel('Imaginary($\omega$)', 'Interpreter','latex')


figure(12);
subplot(1,2,1)
plot(in2, lambdaAngleIn, '.r'); hold on
plot(out2, lambdaAngleOut, '.k');
xlabel('mode', 'Interpreter','latex'); ylabel('$|$angle($\lambda$)$|$ [$0$ $\pi$]', 'Interpreter','latex')
subplot(1,2,2)
plot(f2(in2), lambdaAngleIn, '.r'); hold on
plot(f2(out2), lambdaAngleOut, '.k');
xlabel('Frequency (Hz)', 'Interpreter','latex'); ylabel('$|$angle($\lambda$)$|$ [$0$ $\pi$]', 'Interpreter','latex')






























































