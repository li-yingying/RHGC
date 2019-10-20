% Task: compute initialization z(0), method 2: track  A theta_t + B xi_t
% Let z_{t+1}=\tilde theta_{t+1}^{Index}, \tilde theta_{t+1}=Atheta_t + B xi_t
% Note: this is for GD, if for RHGD, you have to update it online

function z_ini = initial2_tracktildetheta(A, B, theta,xi, N)

[cano_flag, Index, p, n,m,pi] = check_cano(A, B); % Index stores ki


if cano_flag ==0
    disp("Error: A, B not cano")
    z_ini=0;
    return
end

z_ini = zeros(m,N); 

for t=1:N
    % Let z_{t+1}=\tilde theta_{t+1}^{Index}, \tilde theta_{t+1}=Atheta_t + B xi_t
    thetatilde = A*theta(:,t)+ B*xi(:,t);
    z_ini(:,t) = thetatilde(Index); %note theta_t t=1,..., N+1 corresponds to paper's index: theta_0, ..., theta_N
end
