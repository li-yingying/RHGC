% Task: compute initialization z(0), method 1: track previous theta t
% Let z_{t+1}=theta_t^{Index}
% Note: this is for GD, if for RHGD, you have to update it online

function z_ini = initial1_tracktheta(A, B, theta, N)

[cano_flag, Index, p, n,m,pi] = check_cano(A, B); % Index stores ki


if cano_flag ==0
    disp("Error: A, B not cano")
    z_ini=0;
    return
end

z_ini = zeros(m,N); 

for t=1:N
    % let zt= theta_{t-1}^{Index}
    z_ini(:,t) = theta(Index,t); %note theta_t t=1,..., N+1 corresponds to paper's index: theta_0, ..., theta_N
end
