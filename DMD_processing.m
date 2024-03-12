function [X1, X2, Out] = DMD_processing(rawInputSignal, r, matParam, dt)
%DMD_processing - Perform Dynamic Mode Decomposition (DMD) processing on raw input signal.
%
%   [X1, X2, Out] = DMD_processing(rawInputSignal, r, matParam, dt)
%
%   This function performs Dynamic Mode Decomposition (DMD) processing on a raw input signal.
%
%   Dependencies:
%   - HankelMatrix (Function for constructing Hankel matrices)
%   - DMDselective (Function for selective DMD decomposition)
%   - formReconstructTimeDelayed (Function for reconstructing signal from DMD matrix)
%
%   Parameters:
%   - rawInputSignal (matrix): Raw input signal matrix with each row representing a channel.
%   - r (vector): Vector specifying the range of modes to be retained during DMD decomposition.
%   - matParam (struct): Struct containing matrix parameters for Hankel block construction.
%       - NumObsLag (integer): Number of observations to use in constructing Hankel matrix.
%       - Lag (integer): Lag parameter for forming Hankel matrix.
%   - dt (float): Time step between observations in the raw input signal.
%
%   Returns:
%   - X1 (matrix): First set of input data for DMD.
%   - X2 (matrix): Second set of input data for DMD.
%   - Out (struct): Output struct containing DMD-related results.
%       - Phi (matrix): DMD mode matrix.
%       - omega (vector): discrete-time DMD eigenvalues.
%       - lambda (vector): continuous-time DMD eigenvalues.
%       - b (vector): DMD mode amplitudes.
%       - Xdmd (matrix): Reconstructed DMD matrix.
%       - DMDoutputs (matrix): Reconstructed signal from Xdmd.

% Hankel matrix construction
type = "column";
[~, NumObs] = size(rawInputSignal);
n = matParam.NumObsLag;
m = NumObs - n;
disp(['Hankel matrix: Constructing Hankel matrix from signal with ' num2str(NumObs) ' observations...'])
[X1, X2] = HankelMatrix(rawInputSignal, n, m, type);
disp('Done.')
disp(['Hankel matrix is ' num2str(size(X1, 1)) ' rows by ' num2str(size(X1, 2)) ' columns.'])

% DMD decomposition
disp(['DMD: Decomposing signal into ' num2str(r(2)-r(1) +1) ' modes...'])
[Out.Phi, Out.omega, Out.lambda, Out.b, Out.Xdmd] = DMDselective(X1, X2, r, dt);
disp('Done.')

% Signal reconstruction from Xdmd
% Out.DMDoutputs = formReconstructTimeDelayed(Out.Xdmd, NumChannels, NumObs);
% disp('Reconstruction: Done.')
end