% Task: Plot suboptimal MPC, MPC, RHGM under three schemes
% Task: run stable A, and unstable A
% suboptimal MPC: Iter =1, 5 10
% RHGM: initial method 2: track theta_tilde True stepsize
% save workspace

%% convert A, B to canonical form, convert Qt, Rt as well
% to be written 

%% Setup requirement: A, B are canonical. Qt pd, Rt psd

T=30;  % total time
h= 1; % discretization grid size of continuous time
N=T/h; %  number of time slots


% canonical A, B
% A= [0 1
% -1 2];
% B=[0
%     1];

% stable
% A=[0 1
%     -1/6 5/6];
% B=[0
%     1];

A= [0 1
    1 1];
B=[0
    1];
% A=[0 1 0
%     0 0 1
%     2 1 3];
% 
% B=[0;0;1];

eig(A)

% A=[0 1 0
%     3 0.5 1
%     3 0.5 0];
% B=[ 0 0
%     1 0 
%     0 1];  
% 

% A = [0 1 0 0 
%     0 0 1 0
%     1 0 0 0
%     1 0 0 0];
% B=[ 0 0
%     0 0
%     1 0
%     0 1];




[cano_flag, Index, p, n,m,pi] = check_cano(A, B);


if cano_flag ==0
    disp("Error: A, B not cano")
    return
end

mrA=A(Index,:);
x0=zeros(n,1);

Q = zeros(n, n, N+1);  % all you need for ft: ft, ..., f_{t+p-1}
theta = zeros(n,N+1); 


% EMP: R(:,:,i) start at t-1, Q(:,:,i) start at t+i, so for t, it is Q(:,:,i+1)
% R(:,:,i+2)
R= zeros(m,m, N);  % all you need for g: g_{t-1}, ..., g_{t+p-1}
xi = zeros(m, N);


% pick some costs, remember that Q, R should be pd

% EMP: Q_t => Q(:,:,t+1), R_t ==> R(:,:,t+1), so is theta and xi
% compute mu and lf 
mu =1000;
lf =0;
for i=1:N+1
    q = rand(n,1)+1;
    %q=i*ones(n,1);
    if i<N+1
    %q=ones(n,1);
    Q(:,:,i) = diag(q);
    else
        q=(p+1)*ones(n,1);
        Q(:,:,i)=diag(q);
    end
    
    
    Q_eig = eig(Q(:,:,i));
    mu = min(min(mu, Q_eig));
    
    lf = max(max(lf, Q_eig));
    
    theta(:,i)= 10*rand(n,1)*2-5;
    %theta(:,i)=i;
    %theta(:,i)= 1;
    %theta(:,i)= 10*sin(i*2)-5;
end


% compute lg
lg=0;
for i=1:N
    r = rand(m,1)+1;
    %r= (i+4)*ones(m,1);
    %r= ones(m,1);
    R(:,:,i) = diag(r);
    R_eig = eig(R(:,:,i));
    lg = max(max(lg, R_eig));
    
    
    xi(:,i)= 10*rand(m,1)*2-5;
    %xi(:,i)=2*i+10;
    %xi(:,i)=2;
end



W =15;


%% compute Hessian
Hessian = zeros(m*N);
 [b0, ~,~] = gradient_Cz(A, B, Q, R,theta, xi, x0,zeros(m,N),N, Index, mrA, pi,n,m);
 for i =1:m*N
     ei = zeros(m*N,1);
     ei(i)=1;
     eimatrix = reshape(ei,[m,N]);
     [hi, ~,~] = gradient_Cz(A, B, Q, R,theta, xi, x0,eimatrix,N, Index, mrA, pi,n,m);
     Hessian(:,i) =hi-b0;
 end
 
eiglist = eig(Hessian);
L = max(eiglist);
lambda = min(eiglist);

zeta = L/lambda;
zeta


%{
%% estimate mu, L, zeta, 
% use my paper's method to estimate (use better estimation later)

% kappa
ImrA = [eye(m), -mrA];
kappa= norm(ImrA,2)^2;

%kappa =kappa;
% lf = max(max(eig(Qt_cano))+max(eig(PN_cano)));
% lg= Rt_cano;

L = p*lf+(p+1)*lg*kappa; % not true Lip factor, I just use OCO's formula to guess it
 
% 
 
%
lambda =mu; 
% 
% 
 zeta= L/lambda;
% 
% 
 %epsilon = (1+phi)*h;
 %beta = phi^2/(2-phi);
 %gamma = phi^2/((1+phi)*(2-phi));
 %delta = phi^2/(1-phi^2);

%}

%% optimal cost
[Copt, xopt, uopt,  uopt_const,  uopt_gain] = opt_control_LQT(A, B, Q, R,theta, xi, n, m, x0,N);


%
%% RHGD

% stepsize


% stepsize
epsilon =1/L;
beta =0; gamma =0;delta =0;

% solve RHTM


xrhgd = zeros(n,N+1,W);
urhgd = zeros(m,N+1,W);
Crhgd = zeros(W,1);

for wind = 1:W
    xrhgd(:,1,wind)=x0;
    [~, xrhgd(:,2:end,wind),urhgd(:,2:end,wind), Crhgd(wind),gradz_rhgd,K, z_ini_rhgd] = RHGM_ver2_initial1(A, B, wind,  epsilon, beta, delta, gamma, x0, Q, R, theta, xi,N);
end
RegretRHGD = Crhgd-Copt;


%% RHAG

% stepsize

% stepsize
epsilon =1/L;
beta =(sqrt(zeta)-1)/(sqrt(zeta)+1); gamma =beta;delta =0;

% solve RHTM


xrhag = zeros(n,N+1,W);
urhag = zeros(m,N+1,W);
Crhag = zeros(W,1);

for wind = 1:W
    xrhag(:,1,wind)=x0;
    [~, xrhag(:,2:end,wind),urhag(:,2:end,wind), Crhag(wind),gradz_rhag,K, z_ini_rhag] = RHGM_ver2_initial1(A, B, wind, epsilon, beta, delta, gamma, x0, Q, R, theta, xi,N);
end
RegretRHAG = Crhag-Copt;


%% RHTM

% stepsize
 phi= 1-1/sqrt(zeta);

epsilon = (1+phi)/L;
 beta = phi^2/(2-phi);
 gamma = phi^2/((1+phi)*(2-phi));
 delta = phi^2/(1-phi^2);
% solve RHTM




xrhtm = zeros(n,N+1,W);
urhtm = zeros(m,N+1,W);
Crhtm = zeros(W,1);

for wind = 1:W
    xrhtm(:,1,wind)=x0;
    [~, xrhtm(:,2:end,wind),urhtm(:,2:end,wind), Crhtm(wind),gradz_rhtm,K, z_ini_rhtm] = RHGM_ver2_initial1(A, B, wind,  epsilon, beta, delta, gamma, x0, Q, R, theta, xi,N);
end

RegretRHTM = Crhtm-Copt;



%% suboptimal MPC


%== Iter 1
Iter =1;
Csubmpciter1=zeros(W,1);

for wind =1:W
[x,u, Csubmpciter1(wind)] = suboptMPC_morari(A, B, Q, R, theta,xi, N, wind, x0,Iter);
end

RegretSubiter1 = Csubmpciter1-Copt;



%== Iter 5
Iter =5;
Csubmpciter5=zeros(W,1);

for wind =1:W
[x,u, Csubmpciter5(wind)] = suboptMPC_morari(A, B, Q, R, theta,xi, N, wind, x0,Iter);
end

RegretSubiter5 = Csubmpciter5-Copt;


%%

%== Iter 10
Iter =3;
Csubmpciter3=zeros(W,1);

for wind =1:W
[x,u, Csubmpciter3(wind)] = suboptMPC_morari(A, B, Q, R, theta,xi, N, wind, x0,Iter);
end

RegretSubiter3 = Csubmpciter3-Copt;





%% plot

figure;
list =1:W;
plot(1:W, log(RegretRHGD),1:W, log(RegretRHAG), 1:W, log(RegretRHTM),list, log(RegretSubiter1(list)), ...
    list, log(RegretSubiter3(list)), list, log(RegretSubiter5(list)));
legend("RHGD","RHAG","RHTM", "subMPC Iter = 1", "subMPC Iter = 3", "subMPC Iter = 5")
title("log(regret)")

%% plot not first several W
list =3:W;
figure;
plot(list, log(RegretRHGD(list)),list, log(RegretRHAG(list)), list, log(RegretRHTM(list)),list, log(RegretSubiter1(list)), ...
    list, log(RegretSubiter3(list)), list, log(RegretSubiter5(list)));
legend("RHGD initial3","RHAG initial3","RHTM initial3", "subMPC Iter = 1", "subMPC Iter = 3", "subMPC Iter = 5")
title("log(regret)")

















