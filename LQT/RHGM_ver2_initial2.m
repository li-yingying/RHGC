
% task: implement RHGM (triple momentum), for LQT, time-varying Qt, Rt,
% initialization is method 2: track  A theta_t + B xi_t:
% Let z_{t+1}=\tilde theta_{t+1}^{Index}, \tilde theta_{t+1}=Atheta_t + B xi_t
% thetat, xit
% At 2:39pm May 5 2019: this is the correct version, old version RHGM.m is wrong


function [z, x,u, Cost,gradz,K, z_ini] = RHGM_ver2_initial2(A, B, W,  epsilon, beta, delta, gamma,  x0, Q, R, theta, xi,N)

% === input ===
% input: canonical form A, B, p, W, K, epsilon, beta, delta, gamma, z_ini
% for z(0) for all t, 
% Index: \mathcal I: rows of A with *, N horizon
% x0 (initial state) (may not be 0), Qt, Rt, thetat,xit for all t

% === output ===
% z(k) for k=1:K, u(k) for k=1:K, for all t, cost for all k
% z: 1:N, do not have z0



% === procedure: test canonical form (later)

[cano_flag, Index, p, n,m,pi] = check_cano(A, B); % Index stores ki


if cano_flag ==0
    disp("Error: A, B not cano")
    z=0;
    x=0; u=0; Cost=0; gradz=0;K=0;
    return
end


% preparation

mrA = A(Index,:);
K = floor((W-1)/p);
z_ini = zeros(m,N);
y_ini = zeros(m,N);
omega_ini0 = zeros(m,N);
omega_ini1 = zeros(m,N);


    u = zeros(m,N); % compute control at z(K). Note: you can also compute control for all K, but that's for later evaluation.
    x= zeros(n,N); % EMP: x does not include x0

    
% === compute z_ini_1 
%{ 
z_ini_1 is z's value for t<=0, determined by x0, 
 constant for all k
%}

z_ini_1= zeros(m,p); % z_{1-p},..., z_0
for i =1:m
    for j=0: pi(i)-1
        z_ini_1(i, -j-(1-p)+1)= x0(Index(i)-j);
    end
end



% record gradient, for z1,..., zN
gradz = zeros(m,N,K); 


% if K==0 then output initial value
if K==0
    
    for t= 1-W:N-1
        
        % initialize z_ini(:,t+W) if t+W <=N
        if t+W<=N
            thetatilde = A*theta(:,t+W)+ B*xi(:,t+W);
            z_ini(:,t+W) = thetatilde(Index); %note theta_t t=1,..., N+1 corresponds to paper's index: theta_0, ..., theta_N

        end
    
        
        % compute x u
        % EMP: u(:,t+1)=ut, x(:,t)=xt, z(:,t,k) =z_t(k)
        if t==0
            u(:,t+1) = z_ini(:,t+1)-mrA*x0;  % control at t is stored at index t+1
            x(:,t+1)= A*x0+B*u(:,t+1);
        end
        
        if t>=1
            u(:,t+1) = z_ini(:,t+1)-mrA*x(:,t);
            x(:,t+1)=A*x(:,t)+B*u(:,t+1);
        end
    end
    
    z=0; % no update
    Cost = total_cost([x0,x],u,Q,R, theta,xi,N);
    
    
else  % if K>=1
    
    
    
   
    
    
    
    
    % EMP: note: t=1,..., N, so y_t(k)=y(:,t,k), same for omega, and z,
    % different from Q, R
    omega = zeros(m, N,K); % omega_t (k) K iteration, k=1,..., K without initials
    y=zeros(m,N, K);  % note: t=1,..., N, so y_t(k)=y(:,t,k)
    z=zeros(m,N,K);
    
    
    
    for t=1-W:N-1
        
        
        % === initialize z_ini at t+W, also y_ini, omega_ini0, omega_ini1
        
        % initialize z_ini(:,t+W) if t+W <=N
        if t+W<=N
            thetatilde = A*theta(:,t+W)+ B*xi(:,t+W);
            z_ini(:,t+W) = thetatilde(Index); %note theta_t t=1,..., N+1 corresponds to paper's index: theta_0, ..., theta_N
            y_ini(:,t+W) = z_ini(:,t+W);  % initialize omega(0), omega(-1), y(0)
            omega_ini0(:,t+W) = z_ini(:,t+W); % omega(0),
            omega_ini1(:,t+W) = z_ini(:,t+W);    % omega(-1),
        end
        
        
         
         % === To be changed: only update t+W
        
    
    
        
    % == update ==
        
        % t+W-Kt p >=1
        Kt = floor( ( min(t+W, W)-1)/p);
        for j = 1:Kt
            t2 = t+W-j*p; % t2 is global (paper time), t2=1,..., N
            
            if t2<=N  % we also have t2>=1 by def of Kt
                
                % == compute partial gradient partial C partial z_{t2} at y(j-1)
                
                % generate Q_short: t2,..., t2+p-1, if tau >N, Qtau=0
                Q_short = zeros(n,n,p);
                
                for tau= t2: min(t2+p-1, N)  % tau is global time index, tau +1 is Q's time index
                    Q_short(:,:,tau-t2+1)=Q(:,:,tau+1);  % tau-t2+1 is Q_short's time index
                end
                
                
                % generate theta_short
                theta_short = zeros(n,p);
                
                for tau= t2: min(t2+p-1, N)  % tau is global time index
                    theta_short(:,tau-t2+1)=theta(:,tau+1);  % tau-t2+1 is Q_short's time index
                end
                
                % generate R_short: R_{t-1},..., R_{t+p-1}, if tau>N-1, Rtau=0
                R_short = zeros(m,m,p+1);
                
                for tau= t2-1: min(t2+p-1, N-1)  % tau is global time index, tau +1 is R's time index
                    R_short(:,:,tau-(t2-1)+1)=R(:,:,tau+1);  % tau-t2+1 is R_short's time index
                end
                
                % generate xi_short
                xi_short = zeros(m,p+1);
                
                for tau= t2-1: min(t2+p-1, N-1)  % tau is global time index, tau +1 is R's time index
                    xi_short(:,tau-t2+2)=xi(:,tau+1);  % tau-t2+1 is R_short's time index
                end
                
                
                % generate y_short: y_{t-p}, ..., y_{t+p}. If tau>N, let it be 0, cause
                % its value won't affect gradient. If tau<=0, let it be equal to
                % z_ini_1
                
                
                
                y_short= zeros(m, 2*p+1);  % this is the result of all initial values are 0
                if t2-p>=1
                    for tau = t2-p:min(t2+p, N)  % tau is the global time index
                        if j==1 % initial
                            y_short(:,tau-(t2-p)+1) = y_ini(:,tau);  % tau-(t2-p)+1 is the time index in y_short
                        else
                            y_short(:,tau-(t2-p)+1) = y(:,tau, j-1);
                        end
                    end
                end
                
                
                
                if t2-p<=0
                    for tau = t2-p:0
                        y_short(:,tau-(t2-p)+1) = z_ini_1(:,tau-(1-p)+1);  % tau-(1-p)+1 is the time index in z_ini_1
                    end
                    
                    
                    for tau = 1:min(t2+p, N)
                        if j==1 % initial
                            y_short(:,tau-(t2-p)+1) = y_ini(:,tau);
                        else  % j-1 th update
                            y_short(:,tau-(t2-p)+1) = y(:,tau, j-1);
                        end
                    end
                end
                
                
                
                % == compute partial grad
                % Require: y_short's length is (m, 2p+1)
                
               
                if j>1 && t2>1
                    t;
                end
                    
                [pgrad,~,~, pgradf, pgradg, pgradg_1] = grad_eval_LQT(A,p, m,n, Index, y_short,Q_short,R_short,theta_short, xi_short,t2);
                
                if j>1
                    pgradf;
                    pgradg;
                    pgradg_1;
                end
                
                % == pgrad is gradient wrt zt2 at y(j-1)
                gradz(:,t2,j)=pgrad;
                
                % == update by TM
                if j ==1
                    omega(:,t2,j) =(1+beta)* omega_ini0(:,t2)-beta*omega_ini1(:,t2)-epsilon*pgrad;
                    y(:,t2,j) = (1+gamma)*omega(:,t2,j)-gamma*omega_ini0(:,t2);
                    z(:,t2,j) = (1+delta)*omega(:,t2,j)-delta*omega_ini0(:,t2);
                end
                
                
                if j ==2
                    omega(:,t2,j) =(1+beta)* omega(:,t2,j-1)-beta*omega_ini0(:,t2)-epsilon*pgrad;
                    y(:,t2,j) = (1+gamma)*omega(:,t2,j)-gamma*omega(:,t2,j-1);
                    z(:,t2,j) = (1+delta)*omega(:,t2,j)-delta*omega(:,t2,j-1);
                end
                
                
                
                if j >2
                    omega(:,t2,j) =(1+beta)* omega(:,t2,j-1)-beta*omega(:,t2,j-2)-epsilon*pgrad;
                    y(:,t2,j) = (1+gamma)*omega(:,t2,j)-gamma*omega(:,t2,j-1);
                    z(:,t2,j) = (1+delta)*omega(:,t2,j)-delta*omega(:,t2,j-1);
                end
            end
        end
        
        
        % == compute control
        
        
        % EMP: u(:,t+1)=ut, x(:,t)=xt, z(:,t,k) =z_t(k)
        if t==0
            u(:,t+1) = z(:,t+1,K)-mrA*x0;  % control at t is stored at index t+1
            x(:,t+1)= A*x0+B*u(:,t+1);
        end
        
        if t>=1
            u(:,t+1) = z(:,t+1,K)-mrA*x(:,t);
            x(:,t+1)=A*x(:,t)+B*u(:,t+1);
        end
        
        
    end
    Cost = total_cost([x0,x],u,Q,R, theta,xi,N);
end

















