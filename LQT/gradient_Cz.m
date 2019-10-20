% Task: compute gradient of C(z) at z
% full gradient of C(z)

function [gradz, x_fromz, u_fromz] = gradient_Cz(A, B, Q, R,theta, xi, x0, z,N, Index, mrA, pi,n,m)
% A, B canonical form
% pi = [p1,..., pm], 
% Q: t=0: N so is theta
% R: t=0:N-1 so is xi
% z: m by N: z1,..., zN
% x0: initial state
% output: gradz: grad of C(z) at z
% x_fromz: x1,..., xN n by N. $EMP: x_fromz does not have x0
% u_fromz: u0, ..., uN-1, m by N

% EMP: Q1, theta 1 useless



%% compute Mx, hx, Mu, hu such that x=Mx z +hx, u= Mu z + hu
nx = n*N;  % x: x1,..., xN
nu = m*N;  % u: u0,..., uN-1
nz = m*N; % z: z1,..., zN
Mx = zeros(nx, nz);
hx = zeros(nx,1);
Mu = zeros(nu,nz); hu = zeros(nu,1);
ki = Index;

% compute Mx, hx

for t= 1:N
    for i=1:m
        for j=0:pi(i)-1
            if t-j>=1
                Mx(n*(t-1)+ki(i)-j, (t-j-1)*m+i)=1;
            end
            if t-j<=0
                hx(n*(t-1)+ki(i)-j)=x0(ki(i)-j+t);
            end
        end
    end
end


% % use symbolic to test
% x0sym = sym('xi%dt0',[n,1]);
% hxsym = zeros(nx,1); % compute hx symbol
% hxsym=sym(hxsym);
% for t= 1:N
%     for i=1:m
%         for j=0:pi(i)-1
%            
%             if t-j<=0
%                 hxsym(n*(t-1)+ki(i)-j)=x0sym(ki(i)-j+t);
%             end
%         end
%     end
% end
% 
% 
% 
% xbf = sym('xi%dt%d',[n,N],'real');
% zbf = sym('zi%dt%d',[m,N],'real');
% zbfvec= zbf(:);
% xbfvec= xbf(:);
% 
% Mxsym = sym(Mx);
% xfromz = Mxsym*zbfvec+hxsym;
% 


% compute Mu
Iz = eye(nz);

Mu1 = zeros(nu, nx);
for i = 2:N
    Mu1((i-1)*m+1:i*m,(i-2)*n+1:(i-1)*n)=mrA;
end

Mu = Iz -Mu1*Mx;


% compute hu
hu1= zeros(nu,1);
hu1(1:m)= mrA*x0;
hu2= Mu1*hx;
hu=-hu1-hu2;


% 
% % use symbolic to test
% Musym=sym(Mu);
% hu1sym= zeros(nu,1);
% hu1sym=sym(hu1sym);
% mrAsym= sym(mrA);
% hu1sym(1:m)= mrAsym*x0sym;
% 
% Mu1sym=sym(Mu1);
% hu2sym= zeros(nu,1);
% hu2sym=sym(hu2sym);
% 
% husym= zeros(nu,1);
% husym=sym(husym);
% 
% hu2sym= Mu1sym*hxsym;
% husym=-hu1sym-hu2sym;
% 
% ufromz=Musym*zbfvec+husym
% 
%     


%% grad of z


zvec = z(:);
xvec = Mx*zvec + hx; % the corresponding x of z
uvec = Mu*zvec +hu; % the corresponding u of z

x_fromz =reshape(xvec, [n,N]); % n by N x: x1,...,xN $EMP: x u time index different
u_fromz = reshape(uvec, [m,N]); % m by N u: u0,..., uN-1 EMP

gradx= zeros(n*N,1);
gradu=zeros(m*N,1);

for t=2:N+1 % wrt x:1 :N
    gradx((t-2)*n+1:(t-1)*n)= Q(:,:,t)*(x_fromz(:,t-1)-theta(:,t)); % t=2==> Q1(x1-theta1)
end

for t=1:N
    gradu((t-1)*m+1:t*m)= R(:,:,t)*(u_fromz(:,t)-xi(:,t)); %t=1: R0(u0-xi0)
end

gradz= Mx'*gradx+Mu'*gradu;



    






