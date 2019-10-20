% Task: write NAG for H(u), (LQT in terms of u)

function [unext, yvec, uvec,Qf] =NAG_Hu(A, B, Q, R, theta, xi, W, x0, u_ini,Iter)

%{
Input: 
Q t: t+W, but Qt+W is not true, is terminal, so is theta
R: t:t+W-1, so is xi
x0 is xt,
u_ini is m by W, initial state, 
Iter: how many iterations of NAG
%}

%% evaluate Hessian and stepsize

m = size(B,2);

 % == compute Hessian
    Hessian = zeros(m*W);
    [~, b0] = grad_eval_Hu(A, B, Q, R, theta, xi,W,x0, zeros(m,W));
    
 for i =1:m*W
     ei = zeros(m*W,1);
     ei(i)=1;
     eimatrix = reshape(ei,[m,W]);
      [~, hi] = grad_eval_Hu(A, B, Q, R, theta, xi,W,x0,eimatrix);
      Hessian(:,i) =hi-b0;
 end
 
eiglist = eig(Hessian);
LH = max(eiglist);
muH = min(eiglist(find(eiglist>0)));

Qf=LH/muH;


LH;
muH;
Qf;

%% run update

% update K steps
K=Iter;

y_ini = u_ini; % update auxiliary



nu = m*W;

uvec = zeros(nu,K+1);
yvec = zeros(nu, K+1);

uvec(:,1)=u_ini(:);
yvec(:,1) = y_ini(:);

epsilon = 1/LH;



for k=2:K+1
    [~, graduvec] = grad_eval_Hu(A, B, Q, R, theta, xi,W,x0, reshape(yvec(:,k-1),[m,W]));
    uvec(:,k)= yvec(:,k-1)-epsilon* graduvec;
    yvec(:,k)= uvec(:,k) + (sqrt(Qf)-1)/(sqrt(Qf)+1)*(uvec(:,k)  - uvec(:,k-1) );
end


unext =reshape(uvec(:,K+1),[m,W]);














