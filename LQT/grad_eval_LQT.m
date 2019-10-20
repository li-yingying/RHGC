

% goal: compute $\frac{\partial C}{\partial z_t}(z)$

function [pgrad,x1,x,pgradf, pgradg, pgradg_1] = grad_eval_LQT(A,p, m,n, Index, z,Q,R,theta, xi,t)

% input: ft, ..., ft+p-1, g_{t-1},... g_{t+p-1}, z_{t-p},... z_{t+p},
% global t
% output: partial C/partial z_t

%% == input ==
% relevant costs: Q=(Qt, ..., Q_{t+p-1}), R=(R_{t-1},..., R_{t+p-1}), theta=(theta_t, ...,
% theta_{t+p-1}), xi =(xi_{t-1},..., xi_{t+p-1})

% relevant z: z=(z_{t-p}, ..., z_{t+p})
% current time: t

% structure of A, B canonical form: m, n, p, A, Index is \mathcal I={k1,
% ..., km) where ki is A's rows with *

%% == Output ==
% pgrad = $\frac{\partial C}{\partial z_t}(z)$
% auxiliary outputs: x1=x_{t-1}, x={xt, ..., x_{t+p-1}) based on given z


mrA = A(Index, :);




pgrad =zeros(m,1);  % output partial gradient of C wrt zt
pgradf = zeros(m,p); % record partial gradient of all ft, ..., ft+p-1


x = zeros(n,p); % state based on given z of xt, ..., x_{t+p-1}
gradf = zeros(n,p); % full grad df/dx for all t,..., t+p-1


% find Index such that ki-j<=k_{i-1}
    Index_m1 = [0, Index(1:end-1)];
    
    
    
% compute f's grad
for i=1:p
    
    
    
    % == compute x_{t+i-1}  (% possibly slow) == 
    t_ind = t+i-1;  % current t
    
    % what I need from z: z_{t-p+1}, ..., z_t
    t1= t_ind -p+1;  % global index of t_ind's required z
    t2= t_ind;  % global index
    z_ind_1= t1-(t-p)+1;  % translate to z's index
    z_ind_2 = t2-(t-p)+1; 
    
    
    
    
    % compute x_{t+i-1} by z
    xt_ind= zeros(n,1);
    for j=0:p-1
        Index_z= Index -j;  % update x, needed index I-j
        
        % === I changed here, not sure whether it is correct it use to be
        % find(Index_z>=1)
        Index_x = Index_z(find(Index_z> Index_m1));   % only need entries that satisfy ki -j>k_{i-1}
        xt_ind(Index_x)= z(find(Index_z>Index_m1), z_ind_2-j);   % xt^{I-j} = zt_j
    end
    x(:,i)=xt_ind;
    
    
    
    
    
    % == compute full gradient df_{t+i-1}/dx_{t+i-1} ==
    gradf(:,i)= Q(:,:,i)*(xt_ind-theta(:,i));
    
    
    
    
    
    
    % == compute partial grad ==
    
    % partial f_{t+i-1}/partial z_t
    Index_i = Index - i +1;
    
    % nonpositive index
    nonpo = find(Index_i<=Index_m1);
    dz=ones(m,1);
    dz(nonpo)=0;   % if index_i<=0, then partial grad =0
    DZ = diag(dz);   % multiply it
    
    reIndex_i = max(Index_i, 1); % let index be at least 1, to avoid error of calling negative index
    
    
    pgradf(:,i) = DZ* gradf(reIndex_i,i);
    
    
    % == sum up partial grad for f ==
    pgrad = pgrad + pgradf(:,i);
end

pgradf;



%% compute partial gradient of g
    
pgradg_1 = zeros(m,1);  % partial grad of g_{t-1}
% compute g's grad wrt g_{t-1}

% what's zt's index in z? A: p+1
t_indz= t-(t-p)+1;

% == compute xt-1 from z ==
t_ind = t-1; % current index: t-1

% what I need from z: z_{t-p+1}, ..., z_t
    t1= t_ind -p+1;  % global index of t_ind's required z
    t2= t_ind;  % global index
    z_ind_1= t1-(t-p)+1;  % translate to z's index
    z_ind_2 = t2-(t-p)+1; 
    
    % compute x_{t+i-1} by z
    xt_ind= zeros(n,1);
    for j=0:p-1
        Index_z= Index -j;  % update x, needed index I-j
        Index_x = Index_z(find(Index_z> Index_m1));   % only need entries in I-j that is ki-j>k_{i-1}
        xt_ind(Index_x)= z(find(Index_z > Index_m1), z_ind_2-j);   % xt^{I-j} = zt_j
    end
    x1=xt_ind; %x1 is x(t-1)
    
    

% == compute partial grad of g_{t-1} wrt zt
pgradg_1= R(:,:,1)*(z(:,t_indz)-mrA*x1-xi(:,1));
pgrad=pgrad+pgradg_1;
pgradg_1;





% == compute g's grad wrt g_t, ..., g_{t+p-1} ==


pgradg = zeros(m,p); % record partial gradient of all g_{t},... g_{t+p-1},

mid_gradg= zeros(n,p);  % gradient of g_{t} wrt xt by z_{t+1}-mr A x_t -xi_t =u_t
for i=1:p
    % current t
    t_ind = t+i-1;
    
    % current t+1's index in z's index system
    t1= t_ind+1-(t-p)+1;
    
    mid_gradg(:,i) = mrA'*R(:,:,i+1)*(mrA* x(:,i)+xi(:,i+1)-z(:,t1));
    
    
    
    % == evaluate partial grad ==
    % partial g_{t+i-1}/partial z_t
    Index_i = Index - i +1;
    
    
    nonpo = find(Index_i<=Index_m1);
    dz=ones(m,1);
    dz(nonpo)=0;   % if index_i<=0, then partial grad =0
    DZ = diag(dz);   % multiply it
    
    reIndex_i = max(Index_i, 1); % let index be at least 1, to avoid error of calling negative index
    
    
    pgradg(:,i) = DZ* mid_gradg(reIndex_i,i);
    
    
    % == sum up partial grad for f ==
    pgrad = pgrad + pgradg(:,i);
end
pgradg;
t;











