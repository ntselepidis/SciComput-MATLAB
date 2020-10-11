function alpha = nystrom_solve(t, sigma, Lambda, U)
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
    
    % Solve linear system (see eq. (11))
    y = Lambda .* (U' * t);

    z = (Lambda .* (U' * U) + sigma * speye(length(Lambda))) \ y;

    alpha = (1.0 / sigma) * (t - U * z);

end
