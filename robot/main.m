
% ------------------------------------------------------------------------ %
% --------------- Online Control with Limited Prediction ----------------- %
% --------------- Two-wheel robot trejectory tracking -------------------- %
% ------------------------------------------------------------------------ %

clear;
close all;
clc;

%% ========================= Parameter Setting ============================

% ------- load reference trajectory
% 1: circle_1s   2: circle_0.1s  3: square_1s    4: square_0.1s   5:  heart

Type_Trajectory = 5;

[dt_con, T_con, t_con, nT_con, xr, yr, ar]  = load_ref_tra(Type_Trajectory );

% ------- cost
c_x = 1;    c_y = 1;      c_a = 1;   % tracking error
c_v = 15*(dt_con)^2;      c_w = 15*(dt_con)^2;     % velocity effort

% ------- step size
Lx = 0.0015;   Ly = 0.0015;

%% ================== Online Control Algorithm ===========================

% ------- prediction window 
W = 40;
K_inner = floor((W-1)/2);   % inner iteration number

% ------- real trajectory 
x_real = zeros(nT_con,1);   y_real = zeros(nT_con,1);    a_real = zeros(nT_con,1);
v_real = zeros(nT_con,1);   w_real = zeros(nT_con,1);

x_real(1) = xr(1)-2;    y_real(1) = yr(1);      a_real(1) = ar(1); % starting point

% ------- initial for online decision variables
x_ol = zeros(1+W,1);   y_ol = zeros(1+W,1);    a_ol = zeros(1+W,1);

% ------- simulation time
dt_sim = 0.001;
nT_sim = floor(dt_con/dt_sim); 
x_sim = zeros(1e5,1);  y_sim = zeros(1e5,1);

% inner real simulation 
x_temp = zeros(nT_sim,1);  y_temp = zeros(nT_sim,1);    a_temp = zeros(nT_sim,1);

x_last(3:W+1) =  xr(2:W);  y_last(3:W+1) =  yr(2:W);

%% Run Online Algorithm
for t = 1:(nT_con-W)
    
    % 1 - initial
    x_ol(1) = x_real(t);      y_ol(1) = y_real(t);  % fixed 

    % initial with the last iteration results
    x_ol(2:W) = x_last(3:W+1);
    y_ol(2:W) = y_last(3:W+1);
    x_ol(end) = xr(t+W);   % to optimize
    y_ol(end) = yr(t+W);
    
    % 2 - gradient 
    for i = 1:K_inner  
        
        W_temp = 3+W-2*i;
        
        dJdx = Cal_dJdx (x_ol(1:W_temp),y_ol(1:W_temp),c_x,c_y,c_v,c_w,dt_con,xr(t:t+W_temp-1));
        dJdy = Cal_dJdy (x_ol(1:W_temp),y_ol(1:W_temp),c_x,c_y,c_v,c_w,dt_con,yr(t:t+W_temp-1));
        
        % update
        x_ol(2:W_temp-2) = x_ol(2:W_temp-2)- Lx*dJdx;
        y_ol(2:W_temp-2) = y_ol(2:W_temp-2)- Ly*dJdy;
        
    end
    
     x_last = x_ol;     y_last = y_ol; % record online decision
    
    % 3 - calculate control decision
    % tangential velocity
    v_real(t) = 1/dt_con*sqrt((x_ol(2)-x_ol(1))^2  + (y_ol(2)-y_ol(1))^2 );
    
    % angular velocity
    w_real(t) = cal_ang_vel( x_ol(1:2), y_ol(1:2), a_real(t),dt_con);
        
    % 4 - real dynamics
    x_temp(1) = x_real(t);    y_temp(1) = y_real(t);   a_temp(1) = a_real(t);
    for tt = 1:(nT_sim-1)
        x_temp(tt+1) =  x_temp(tt) + dt_sim*cos(a_temp(tt))* v_real(t);
        y_temp(tt+1) =  y_temp(tt) + dt_sim*sin(a_temp(tt))* v_real(t);
        a_temp(tt+1) =  a_temp(tt) + dt_sim*w_real(t);
    end
    
    % update
    x_real(t+1) = x_temp(end); 
    y_real(t+1) = y_temp(end); 
    a_real(t+1) = a_temp(end); 
    
    x_sim((1+(t-1)*nT_sim):(t*nT_sim)) = x_temp;
    y_sim((1+(t-1)*nT_sim):(t*nT_sim)) = y_temp;

%      x_real(t+1) = x_ol(2); 
%      y_real(t+1) = y_ol(2); 

end

%% Results Display

% dynamics
fig=figure(1); xlabel('X');  ylabel('Y');  axis([-20,20,-20,15]);
for j= 1:floor(0.58*(nT_con))
    pause(0.01);
    clf(fig); xlabel('X');  ylabel('Y');  axis([-20,20,-20,15]);
    hold on; plot(xr(j),yr(j),'r.','markersize',45);
    hold on;  plot(x_real(j),y_real(j),'b.','markersize',45);

    hold on;  plot(xr(1:j),yr(1:j),'r--','LineWidth',1.8); 
    hold on;  plot(x_real(1:j),y_real(1:j),'b-','LineWidth',1.5); 
end


% final result
num_plot_step = floor(0.56*(nT_con));
figure();
plot(xr(1:num_plot_step),yr(1:num_plot_step),'--r','LineWidth',1.5); 
hold on; 
plot(x_sim(1:num_plot_step*nT_sim),y_sim(1:num_plot_step*nT_sim),'-b','LineWidth',1.2); 
xlabel('X'); ylabel('Y');  axis([-20,20,-20,15]);  legend('reference','real');








