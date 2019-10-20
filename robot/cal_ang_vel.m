function w = cal_ang_vel( x, y, a,dt)
   

    angle_next = atan( (y(2)-y(1))/(x(2)-x(1)));  
    sign_x = sign(x(2)-x(1));  sign_y = sign(y(2)-y(1));
    
    if sign_x >=0
        if sign_y >=0  % 1 area
%             angle_next = angle_next;
        else           % 4 area
            angle_next = 2*pi + angle_next;
        end
    else
        if sign_y >=0   % 2 area
            angle_next = pi + angle_next;
        else            % 3 area
            angle_next = pi + angle_next;
        end
    end
    
    angle_now = mod(abs(a),2*pi); % tramsform to (0,2*pi)
    if a < 0
       angle_now = 2*pi - angle_now;
    end
    
    angle_diff = angle_next - angle_now;
    
    w = 1/dt*( angle_diff );
   




end

