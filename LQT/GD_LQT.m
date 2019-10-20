% Task: LQT, do GD for C(z) return iteration z

function [zvec, zmatrix, gradz, xz, uz, Cost] = GD_LQT(A, B, Q, R, theta, xi,  x0, N,z_ini, epsilon, K)
% input: A, B, cano, so are costs
% x0 initial state, n by 1
% z_ini initial z m by N
% epsilon is stepsize
%K is iteration number
% output: zvec zmatrix: z(k)k=0:K has initial
% gradz: grad(k) k=0:K
% xz, uz(k): k=0:K Note: xz: n by N+1, has x0, 
% Cost(k): k=0:K

[cano_flag, Index, p, n,m,pi] = check_cano(A, B);

if cano_flag ==0
    disp("Error: A, B not cano")
    z=0; 
    gradz=0;
    return
end

mrA=A(Index,:);


zvec= zeros(m*N,K+1); 
zmatrix = zeros(m,N, K+1);
gradz= zeros(m*N,K+1);
zmatrix(:,:,1)=z_ini;
zvec(:,1)=z_ini(:);

Cost = zeros(K+1,1);

xz= zeros(n,N+1, K+1);
uz= zeros(m,N,K+1);
for k =1:K
    xz(:,1,k)=x0;
    
    [gradz(:,k), xz(:,2:end,k), uz(:,:,k)] = gradient_Cz(A, B, Q, R,theta, xi, x0, zmatrix(:,:,k),N, Index, mrA, pi,n,m);
    zvec(:,k+1)= zvec(:,k)-epsilon*gradz(:,k);
    zmatrix(:,:,k+1)= reshape(zvec(:,k+1),[m,N]);
    Cost(k) = total_cost(xz(:,:,k),uz(:,:,k),Q,R, theta,xi,N);
end
k = K+1;
xz(:,1,k)=x0;
    [gradz(:,k), xz(:,2:end,k), uz(:,:,k)] = gradient_Cz(A, B, Q, R,theta, xi, x0, zmatrix(:,:,k),N, Index, mrA, pi,n,m);

    Cost(k) = total_cost(xz(:,:,k),uz(:,:,k),Q,R, theta,xi,N);



