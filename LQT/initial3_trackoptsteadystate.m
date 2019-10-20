% Task: compute initialization z(0), method 3: track  optimal steady state
% Let z_{t+1}= x_{t+1}^{Index}(e), x(e)=Ax(e)+B u(e), 
% Note: this is for GD, if for RHGD, you have to update it online

function [z_ini,xeN,ueN] = initial3_trackoptsteadystate(A, B, Q, R, theta,xi, N)

[cano_flag, Index, p, n,m,pi] = check_cano(A, B); % Index stores ki


if cano_flag ==0
    disp("Error: A, B not cano")
    z_ini=0;
    return
end

z_ini = zeros(m,N); 
xeN= zeros(n,N);
ueN = zeros(m,N);

for t=1:N
    % compute optimal steady state
    
    
    cvx_begin quiet
variable xe(n)
variable ue(m)

minimize((xe-theta(:,t))'*Q(:,:,t)*(xe-theta(:,t))+(ue-xi(:,t))'*R(:,:,t)*(ue - xi(:,t)))

subject to 
xe == A*xe + B*ue;

cvx_end

xeN(:,t)=xe;
ueN(:,t)=ue;


    
    z_ini(:,t) = xe(Index); %note theta_t t=1,..., N+1 corresponds to paper's index: theta_0, ..., theta_N
end
