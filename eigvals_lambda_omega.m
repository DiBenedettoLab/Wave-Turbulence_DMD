function [lambdaComponents, lambdaAngle, circleArc, omegaComponents] = eigvals_lambda_omega(lambda, omega)
% EIGVALS_LAMBDA_OMEGA Extracts components and angles from eigenvalues and omegas
%   [lambdaComponents, lambdaAngle, circleArc, omegaComponents] = eigvals_lambda_omega(lambda, omega)
%   computes the real and imaginary components of eigenvalues and omegas,
%   calculates the absolute angle of eigenvalues, and generates a unit
%   circle arc for better visualization.

% Input:
%   - lambda: discrete -time eigenvalues
%   - omega: continuous-time eigenvalues for predictions not in dt

% Output:
%   - lambdaComponents: Real (row 1) and imaginary (row 2) components of eigenvalues
%   - lambdaAngle: Absolute angle of eigenvalues
%   - circleArc: Coordinates of a unit circle arc for visualization
%   - omegaComponents: Real and imaginary components of omega values

% Initialize output matrices
lambdaComponents = nan(2, size(lambda,1));
lambdaAngle = abs(angle(lambda));
circleArc = nan(2, 500); % 500 points for the unit circle arc
omegaComponents = nan(2, size(omega,1));

% Extract real and imaginary components of eigenvalues
lambdaComponents(1,:) = real(lambda)';
lambdaComponents(2,:) = imag(lambda)';

% Calculate xy projection angle of eigenvalues
xyProjTheta = atan(lambdaComponents(2,:) ./ lambdaComponents(1,:));

% Unit circle
aEps = 0; % epsilon for better visualization of edge modes
a1 = max(xyProjTheta) + aEps;  % max angle
t = linspace(-a1, a1, 500);
x0 = 0; y0 = 0; % center
r = 1; % radius
circleArc(1,:) = x0 + r*cos(t);
circleArc(2,:) = y0 + r*sin(t);

% Extract real and imaginary components of omega values
omegaComponents(1,:) = real(omega)';
omegaComponents(2,:) = imag(omega)';

end
