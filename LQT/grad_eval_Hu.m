%{
Task: compute gradient of H(u) for LQT: compute gradient of H(u) wrt u
H(u) is LQT, translate all x to u

Note: Morari's suboptimal MPC can be applied to linear MPC, not-quadratic
cost.

%}
function [gradu, graduvec] = grad_eval_Hu(A, B, Q, R, theta, xi,W,x0, u)
%{
 Input: Q: Qt,..., Q_{t+W-1}, Pterm, so is theta
           R: Rt, ..., R_{t+W-1}
           u: m by W, where do you want to evaluate
           x0: initial state
%}


n = length(A);
m= size(B,2);

%% compute Mx, hx, M, h
% x= M u
nx = W*n; 
nu = W*m; % u0,..., u_{W-1}
uvec = u(:); 

Mx = zeros(nx, nu);

for k =1:W
    for j=1:k
        Mx( (k-1)*n+1:k*n, (j-1)*m+1:j*m) = (A^(k-j)) * B;
    end
end


hx = zeros(nx,1);
for k=1:W
    hx((k-1)*n+1:k*n) = A^k*x0;
end


xvec = Mx*uvec + hx;  % corresponding xvec of u

x = reshape(xvec, [n,W]); % x in (n,N) form


%% compute J's grad wrt x, wrt u

Fx = zeros(n, W); % grad of J wrt x
Gu = zeros(m,W); % grad of J wrt u

for t=1:W
    Fx(:,t)= Q(:,:,t+1)*(x(:,t)-theta(:,t+1));
    Gu(:,t)= R(:,:,t)*(u(:,t)-xi(:,t));
end
Fxvec = Fx(:);
Guvec = Gu(:);

%% compute H(u)'s grad

graduvec = Guvec + Mx'*Fxvec;

gradu= reshape(graduvec, [m,W]);



























