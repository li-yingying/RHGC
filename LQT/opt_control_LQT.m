% Task: write optimal cost function for LQT with xi, and optimal controller

function [Copt, xopt, uopt,  uopt_const,  uopt_gain] = opt_control_LQT(A, B, Q, R,theta, xi, n, m, x0,N)





%% optimal total cost

% Task: output optimal total cost, optimal controller, optimal trajectory
% of x and u


% input: Q: t=0: N, R t=0:N-1, theta t=0:N, xi t=0:N-1, N x0, A B, n,m
% x0 is column vector (EMP, different from old code)
% note: A, B does not need to be canonical form

%% compute Pt, Mt, betat, zt
 Pt = zeros(n,n,N+1);
 Mt = zeros(n,n,N);
 betat = zeros(N+1,n);  % EMP: betat row is beta
 ztt=zeros(N+1,n); %EMP: ztt row is zt
 
 thetanew=theta'; % thetanew is theta transpose
 
 Pt(:,:,N+1)= Q(:,:,N+1);
 betat(N+1,:)=thetanew(N+1,:);
 ztt(N+1,:)= betat(N+1,:)*Pt(:,:,N+1);
 for t=N:-1:1
     Mt(:,:,t)= Pt(:,:,t+1)-Pt(:,:,t+1)*B*inv(R(:,:,t)+B'*Pt(:,:,t+1)*B)*B'*Pt(:,:,t+1);
     
     Pt(:,:,t)= Q(:,:,t) + A'*Mt(:,:,t)*A;
     
     beta=inv(Q(:,:,t)+A'*Mt(:,:,t)*A)*(Q(:,:,t) *thetanew(t,:)'-A'*Mt(:,:,t)*(B*xi(:,t)-betat(t+1,:)'));
     betat(t,:)=beta';
     ztt(t,:)= betat(t,:)*Pt(:,:,t);
 end


 
 %% simulate x* and u*, save optimal controller
 xopt = zeros(n,N+1);
 uopt=zeros(m,N);
 xopt(:,1)=x0;  % x0 is row vector 
 
 uopt_const = zeros(m, N);  % constant in the optimal controller, determined by xi, theta
 uopt_gain = zeros(m,n, N); % Gain matrix in front of x in the opt controller
 
 for t= 1:N
     Ktx = inv(R(:,:,t)+B'*Pt(:,:,t+1)*B)*B'*Pt(:,:,t+1)*A;
     uopt_gain(:,:,t)= Ktx;
     
     
     Ktz = inv(R(:,:,t)+B'*Pt(:,:,t+1)*B)*B';
     Ktxi = inv(R(:,:,t)+B'*Pt(:,:,t+1)*B)*R(:,:,t);
     
     uopt_const(:,t) = Ktxi*xi(:,t) + Ktz*ztt(t+1,:)'; % EMP: zt is row vector in ztt
     
     uopt(:,t)= -Ktx*xopt(:,t)+uopt_const(:,t);
     xopt(:,t+1)=A*xopt(:,t)+B*uopt(:,t);
 end
 
 
 %% compute optimal cost
 Copt =  total_cost(xopt,uopt,Q,R, theta,xi,N);
 
 
 











