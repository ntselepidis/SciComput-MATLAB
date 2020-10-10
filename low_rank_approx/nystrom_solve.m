function alpha = nystrom_solve(X, t, m, sigma, k, compute_error)
%NYSTROM_SOLVE Solves (K + sigma I) a = t.
%
%   Reference Paper:
%
%   Williams, Christopher KI, and Matthias Seeger.
%   "Using the Nyström method to speed up kernel machines."
%   Advances in neural information processing systems. 2001.
%
%   https://papers.nips.cc/paper/1866-using-the-nystrom-method-to-speed-up-kernel-machines.pdf
%
    if nargin < 6
        compute_error = false;
    end

    n = size(X, 1);

    % Randomly shuffle X
    Xp = X(randperm(n), :);

    % Compute reduced kernel matrix
    Knm = k(Xp, Xp(1:m,:));

    % Compute eigen-decomposition (see eq. (7))
    [U_m, Lambda_m] = eig(Knm(1:m, :), 'vector');

    % Compute approximate eigenvalues (see eq. (8))
    Lambda = (n / m) * Lambda_m;

    % Compute approximate eigenvectors (see eq. (9))
    U = sqrt(m / n) * Knm * (U_m ./ Lambda_m');

    % Solve linear system (see eq. (11))
    y = Lambda .* (U' * t);
    z = (Lambda .* (U' * U) + sigma * speye(m, m)) \ y;
    alpha = (1.0 / sigma) * (t - U * z);

    if compute_error
        % Compute full kernel matrix
        Knn = k(Xp, Xp);
        alpha_true = (Knn + sigma * speye(n,n)) \ t;
        disp(norm(alpha - alpha_true) / norm(alpha_true));
    end
end
