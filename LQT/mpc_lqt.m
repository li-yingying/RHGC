% Task: write MPC for LQT problem
% output control, state, cost

function [xmpc, umpc, Cmpc] = mpc_lqt(A, B, Q, R, theta, xi, n,m, x0, N, W)

xmpc = zeros(n,N+1);
umpc=zeros(m,N);
xmpc(:,1)=x0;


% at each t, solve t:t+W-1, 
% currently, do not consider terminal cost
for t=1:N-W+1 % until t=N-W+1:N
    
    Qmpc =zeros(n,n,W+1);
    thetampc = zeros(n,W+1);
    
    Qmpc(:,:,1:W) = Q(:,:, t:t+W-1);
    Rmpc = R(:,:,t:t+W-1);
    thetampc(:,1:W) = theta(:,t:t+W-1);
    ximpc = xi(:,t:t+W-1);
    if t==N-W+1
        Qmpc(:,:,end)= Q(:,:,N+1); % if t=N-W+1, terminal cost is Q(N+1)=Q_N
        thetampc(:,end)=theta(:,N+1);
    else
        PN=zeros(n,n); % if t<N-W, terminal cost is 0
    end
    
    x0 =xmpc(:,t);
    
    [~, xmpc1, umpc1,  ~,  ~] = ...
        opt_control_LQT(A, B, Qmpc, Rmpc,thetampc, ximpc, n, m, x0,W);
    if t<N-W+1
        umpc(:,t)=umpc1(:,1);
        xmpc(:,t+1)=xmpc1(:,2);
    else
        umpc(:,t:end)=umpc1;
        xmpc(:,t+1:end)=xmpc1(:,2:end);
    end
end




 Cmpc =  total_cost(xmpc,umpc,Q,R, theta,xi,N);


    
    