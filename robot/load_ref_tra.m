function [dt_con, T_con, t_con, nT_con, xr, yr, ar] = load_ref_tra( index )

switch index
    
    % 1: circle_1s 
    case 1
         dt_con = 1; % second
         T_con = 30*pi;  
         t_con = 0:dt_con:T_con;
         nT_con = length(t_con);
         
         radius = 10;      
         wr = 0.1;      vr = radius*wr;
         xr = radius*cos(wr*t_con)';      
         yr = radius*sin(wr*t_con)';      
         ar = wr*t_con' + pi/2;

         disp('load fast circle successfully!');
    case 2
        
         dt_con = 0.1; % second
         T_con = 30*pi;  
         t_con = 0:dt_con:T_con;
         nT_con = length(t_con);
         
         radius = 10;      
         wr = 0.1;      vr = radius*wr;
         xr = radius*cos(wr*t_con)';      
         yr = radius*sin(wr*t_con)';      
         ar = wr*t_con' + pi/2;
          disp('load slow circle successfully!');
    case 3
        dt_con = 1; % second
        T_con = 100; %20*pi+10;     
        t_con = 0:dt_con:T_con;
        nT_con = length(t_con);
        xr = zeros(nT_con,1);  yr = zeros(nT_con,1);    ar = zeros(nT_con,1);
        xr(1:21) = 0:1:20;  
        xr(22:41) = 20*ones(20,1);    yr(22:41) = 1:1:20;       ar(22:41) = pi/2*ones(20,1);
        xr(42:61) = 20 -(1:1:20);  yr(42:61) = 20*ones(20,1);   ar(42:61) = pi*ones(20,1);
        yr(62:81) = 20 -(1:1:20); ar(62:81) = 3/2*pi*ones(20,1);
        xr(82:101) = 1:1:20;  
         disp('load fast square successfully!');
    case 4
        dt_con = 0.1; % second
        T_con = 100; %20*pi+10;     
        t_con = 0:dt_con:T_con;
        nT_con = length(t_con);
        xr = zeros(nT_con,1);  yr = zeros(nT_con,1);    ar = zeros(nT_con,1);
        xr(1:201) = 0:0.1:20;  
        xr(202:401) = 20*ones(200,1);    
        yr(202:401) = 0.1:0.1:20;       
        ar(202:401) = pi/2*ones(200,1);
        xr(402:601) = 20 -(0.1:0.1:20);  
        yr(402:601) = 20*ones(200,1);   
        ar(402:601) = pi*ones(200,1);
        yr(602:801) = 20 -(0.1:0.1:20);
        ar(602:801) = 3/2*pi*ones(200,1);
        xr(802:1001) = 0.1:0.1:20;  
        disp('load slow square successfully!');
         
    case 5
        dt_con = 0.1;
        T_con = 60;
        t_con = 0:dt_con:T_con;
        nT_con = length(t_con);
        t = linspace(-6,6,601);
        xr = 16*(sin(t)).^3;
        yr = 13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t);    
        ar = zeros(nT_con,1);
        for i = 1:nT_con-1
            ar(i) = atan( (yr(i+1)-yr(i))/(xr(i+1)-xr(i)) );
        end
        disp('load heart successfully!');
    otherwise
        disp('please select a correct type number');
end


end

