function [P2, f2, P1, f1] = DMD_spectrum(b ,lambda, dt)
%% DMD_spectrum
% This function computes the DMD spectrum based on the provided DMD results.

% Input:
%   - b: DMD amplitudes.
%   - lambda: DMD eigenvalues.
%   - dt: Time step.

% Output:
%   - P2: Power spectrum for all frequency modes.
%   - f2: All frequency modes.
%   - P1: Power spectrum for the positive frequencies in f2.
%   - f1: Positive frequencies in f2.

% 2-sided spectrum
P2 = abs(b).^2; %(diag(Phi'*Phi)); %dot(Phi, Phi); % norm(Phi)^2 DMD spectrum (check if need to only take the real or imag part?) || phi ||^2
f2 = imag(log(lambda))./(2*pi*dt); %frequencies
%X1
% 1-sided
fidx = find(f2 >= 0);
f1 = f2(fidx);
P1 = 2*P2(fidx);
end