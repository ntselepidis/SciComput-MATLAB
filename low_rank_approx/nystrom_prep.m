function [Lambda, U] = nystrom_prep(X, kernel, m)
%NYSTROM_PREP Implements the preprocessing step of Nystrom approximation method.
%
%   Reference Paper:
%
%   Williams, Christopher KI, and Matthias Seeger.
%   "Using the Nyström method to speed up kernel machines."
%   Advances in neural information processing systems. 2001.
%
%   https://papers.nips.cc/paper/1866-using-the-nystrom-method-to-speed-up-kernel-machines.pdf
%
    
    n = size(X, 1);

    % Compute reduced kernel matrix
    Knm = kernel(X, X(1:m,:));

    % Compute eigen-decomposition (see eq. (7))
    [U_m, Lambda_m] = eig(Knm(1:m, :), 'vector');

    % Compute approximate eigenvalues (see eq. (8))
    Lambda = (n / m) * Lambda_m;

    % Compute approximate eigenvectors (see eq. (9))
    U = sqrt(m / n) * Knm * (U_m ./ Lambda_m');

end
