function [time_dynamics, Xdmd, XdmdMode] = DMD_reconstruct(X,Phi,omega, b, dt)
% DMD_RECONSTRUCT Reconstructs signals using Dynamic Mode Decomposition (DMD) results.
%   [time_dynamics, XdmdMode, Xdmd] = DMD_reconstruct(X, Phi, omega, b, dt) reconstructs
%   signals based on the given DMD results.

% Input:
%   - X: Hankel matrix.
%   - Phi: DMD modes.
%   - omega: DMD frequencies.
%   - b: DMD amplitudes.
%   - dt: Time step.

% Output:
%   - time_dynamics: Time dynamics of the reconstructed signal.
%   - Xdmd: Reconstructed signal.
%   - XdmdMode: Reconstructed signal of each DMD mode.
% code adapted from Keisuke Fujii and Steve Brunton

% x1 = X(:, 1);
% b = Phi \ x1;

r = length(omega);
t = (0:size(X,2))*dt;% time vector
cols = size(X,2) ;

textprogressbar('Reconstructing signal Hankel block from DMD: ')
textprogressbar(0)

if nargout == 3
    % reconstruct mode by mode (slow and requires a lot of memory)
    time_dynamics = zeros(cols,r);
    XdmdMode = zeros(size(X,1), cols, r);
    for iter = 1:cols
        for mode = 1:r
            time_dynamics(iter,mode) = b(mode) * ( exp( omega(mode) * t(iter) ) );
            XdmdMode(:,iter,mode) = real( Phi(:,mode) * ( b(mode) .* exp( omega(mode) * t(iter) ) ) );
        end
        textprogressbar(iter/(cols)*100)
    end
    % Add all modes 
    disp('Adding all modes...')
    Xdmd = sum(XdmdMode,3);
    textprogressbar('Done.')

elseif nargout == 2
    % reconstruct all at once
    time_dynamics = zeros(r, cols);
    % XdmdMode = zeros(size(X,1), cols);
    for iter = 1:cols
        time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
        textprogressbar(iter/(cols)*100)
    end
    Xdmd=Phi*time_dynamics;

else
    error('Usage: [time_dynamics, Xdmd] = DMD_reconstruct(X,Phi,omega,dt) or [time_dynamics, XdmdMode, Xdmd] = DMD_reconstruct(X,Phi,omega,dt)')
end

textprogressbar(100)
fprintf('\n')

end