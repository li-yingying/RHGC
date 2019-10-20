%{
Task: write suboptimal MPC based on Morari 2009 CDC
Nesterov GD, update based on u
only for LQT currently
warm start, let uN=xi_N (same warm start in Morari,who let uN=0 when xi_N=0
%}

function [x,u, cost_submpc,Qf] = suboptMPC_morari(A, B, Q, R, theta,xi, N, W, x0,Iter)

%{

Input: require W<=N, otherwise my code is wrong
%}



%% suboptimal MPC
Qf = zeros(N,1);
n = length(A);
m = size(B,2);
% for simplicity, let terminal cost =0

u = zeros(m,N);  % store u1,..., uN
x = zeros(n,N+1);  %x1,...,xN+1
x(:,1)=x0; % initial
 
ucur = zeros(m,W,N);  % current 
unext = zeros(m,W,N);  % after update


for t=1:N  % time start at 1, unlike paper who start at 0
   
    
    % == prepare cost needed
    % if shorter than W, then fill with 0
    
    % == Let terminal cost =0
    Q_short = zeros(n,n,W+1);   theta_short = zeros(n,W+1);
    if t+W-1 <= N+1
        Q_short(:,:,1:W)= Q(:,:,t:t+W-1);
        theta_short(:,1:W)= theta(:,t:t+W-1);
    else
        Q_short(:,:,1:(N+1-t+1)) = Q(:,:, t:N+1);
        theta_short(:,1:(N+1-t+1))= theta(:,t:N+1);
    end
   
    
    % == even if terminal cost nonzero, this won't change
    R_short = zeros(m,m,W);
    xi_short = zeros(m,W);
    if t+W-1 <=N
        R_short(:,:,1:W)= R(:,:,t:t+W-1);
        xi_short(:,1:W)= xi(:,t:t+W-1);
    else
        R_short(:,:,1:(N-t+1)) = R(:,:, t:N);
        xi_short(:,1:(N-t+1))= xi(:,t:N);
    end
    
    
    % == prepare ucur 
    if t==1
        
        
        % initialize u and y by xi
        
        
        ucur(:,:,t) = xi(:,1:W);  % let u_t = xi_t
       
            
        
    else
        
        
        
        if t+W-1<=N  % only update u_t+W-1 if t+W-1<=N, o.w. it is 0
            ucur(:,W,t)=xi(:,t+W-1);
            
        end
        
        ucur(:,1:W-1,t) = unext(:,2:end,t-1); % warm start
       
        
    end
    % == NAG
    [unext(:,:,t), yvec, uvec,Qf(t)] =NAG_Hu(A, B, Q_short, R_short, theta_short, xi_short, W, x(:,t), ucur(:,:,t),Iter);
    
    
    %[Copt, xopt, uopt,  uopt_const,  uopt_gain] = opt_control_LQT(A, B, Q_short, R_short,theta_short, xi_short, n, m, x(:,t),W);

    u(:,t)=unext(:,1,t);
    x(:,t+1)=A*x(:,t)+B*u(:,t);
end
    
    
    cost_submpc = total_cost(x,u,Q,R, theta,xi,N);







