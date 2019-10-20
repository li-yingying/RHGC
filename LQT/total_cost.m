%% compute J(x,u)
function cost = total_cost(x,u,Q,R, theta,xi,N)
% Qt t=0:N, Rt: t=0:N-1, theta_t: t=0:N, xi_t : t=0:N-1
% x0 initial, x_t: t=0:N, u_t:t=0:N-1
% Require on input: x(n,N+1) all column vector, x0  in x
% u(n,N) column vector, theta(n,N+1), 

% EMP: ut and Rt has same index, but Qt start at t=0, xt starts at t=1

t=0;
cost =0.5*(x(:,1)-theta(:,1))'*Q(:,:,1)*(x(:,1)-theta(:,1));
% x= (x1, x2,..., xN)
for t=1:N
    cost = cost + 0.5*(x(:,t+1)-theta(:,t+1))'*Q(:,:,t+1)*(x(:,t+1)-theta(:,t+1))+0.5* (u(:,t)-xi(:,t))'*R(:,:,t)* (u(:,t)-xi(:,t));
    % note, this is f_t + g_{t-1}
end

end
