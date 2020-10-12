function x = pipelined_bicgstab(A, b, tol, Nmax, M)
%PIPELINED_BICGSTAB Pipelined BiConjugate Gradients Stabilized Method.
%
%   Reference Paper:
%
%   Cools, S. and Vanroose, W., 2017.
%   The communication-hiding pipelined BiCGstab method for the
%   parallel solution of large unsymmetric linear systems.
%   Parallel Computing, 65, pp.1-20.
%
%   https://www.sciencedirect.com/science/article/pii/S0167819117300406
%

n=length(b);

phat = zeros(n,1);
shat = zeros(n,1);
s = zeros(n,1);
z = zeros(n,1);
zhat = zeros(n,1);
v = zeros(n,1);
omega = 0;

tolb = tol*norm(b);

x = zeros(n,1);
r0=b-feval(A,x); r=r0; rhat=feval(M,r); w=feval(A,rhat); what=feval(M,w);
t=feval(A,what); alpha=(r'*r)/(r'*w); beta=0;
for i=1:Nmax
    r_old = r;
    phat = rhat + beta*(phat - omega*shat);
    s    = w    + beta*(s    - omega*z);
    shat = what + beta*(shat - omega*zhat);
    z    = t    + beta*(z    - omega*v);
    q    = r    - alpha*s;
    qhat = rhat - alpha*shat;
    y    = w    - alpha*z;
    qdoty = q'*y;
    ydoty = y'*y;
    zhat = feval(M,z);
    v = feval(A,zhat);
    omega = qdoty/ydoty;
    x = x + alpha*phat + omega*qhat;
    r = q    - omega*y;
    
%     if (norm(r) < tolb)
%         disp(['pipelined bicgstab needed ' ...
%             num2str(i-1) '.5 iterations'])
%         break;
%     end
    
    rhat = qhat - omega*(what - alpha*zhat);
    w    = y    - omega*(t    - alpha*v);
    r0dotr_old = r0'*r_old;
    r0dotr = r0'*r;
    r0dotw = r0'*w;
    r0dots = r0'*s;
    r0dotz = r0'*z;
    what = feval(M,w);
    t = feval(A,what);
    
    if (norm(r) < tolb)
        disp(['pipelined bicgstab converged at iteration ' num2str(i) ...
              ' to a solution with relative residual ' ...
            num2str( norm(r)/norm(b) ) '.'])
        break;
    end
    
    beta = (alpha/omega)*(r0dotr/r0dotr_old);
    alpha = r0dotr/(r0dotw + beta*r0dots - beta*omega*r0dotz);
end

end
