function [t,u,v,r,x,y,psi,U] = Mariner_ZigZag(mariner,x,ui,t_final,t_rudderexecute,h,maneuver)
% ZIGZAG      [t,u,v,r,x,y,psi,U] = zigzag(ship,x,ui,t_final,t_rudderexecute,h,maneuver)
%             performs the zig-zag maneuver, see ExZigZag.m
%
% Inputs :
% 'ship'          = ship model. Compatible with the models under .../gnc/VesselModels/
% x               = initial state vector for ship model
% ui              = [delta,:] where delta=0 and the other values are non-zero if any
% t_final         = final simulation time
% t_rudderexecute = time control input is activated
% h               = sampling time
% maneuver        = [rudder angle, heading angle]. Default 20-20 deg that is: maneuver = [20, 20] 
%                    rudder is changed to maneuver(1) when heading angle is larger than maneuver(2)
%
% Outputs :
% t               = time vector
% u,v,r,x,y,psi,U = time series
%
% Author:    Thor I. Fossen
% Date:      22th July 2001
% Revisions: 15th July 2002, switching logic has been modified to handle arbitrarily maneuvers

if nargin>7 | nargin<6, error('number of inputs must be 6 or 7'); end
if t_final<t_rudderexecute, error('t_final must be larger than t_rudderexecute'); end
if nargin==6, maneuver = [10,10]; end

N = round(t_final/h);               % number of samples
xout = zeros(N+1,9);                % memory allocation

disp('Simulating...')

u_ship=ui;

for i=1:N+1,
    time = (i-1)*h;
    
    psi = x(6)*180/pi;
    r   = x(3);
    
    if round(time)==t_rudderexecute, 
        u_ship(1)=-maneuver(1)*pi/180; 
    end
    
    if round(time) > t_rudderexecute, 
        if (psi<=-maneuver(2) & r<0),
            u_ship(1) = maneuver(1)*pi/180; 
        elseif (psi>=maneuver(2) & r>0),
            u_ship(1) = -maneuver(1)*pi/180;            
        end   
    end
 
    [xdot,U] = feval(mariner,x,u_ship);       % ship model
    
    x = euler2(xdot,x,h);                     % Euler integration
    
    xout(i,:) = [time,x(1:6)',U,u_ship(1)];  
    
    
end

% time-series
t     = xout(:,1);
u     = xout(:,2);
%u1=u.*U;
v     = xout(:,3);
%v1=v.*U;
r     = xout(:,4); 
%r1=(r.*U)/(160.93);
x     = xout(:,5);
y     = xout(:,6);
psi   = xout(:,7);
U     = xout(:,8);
delta_c = xout(:,9);



% plots
figure(1)
plot(x,y,'linewidth',2),grid,axis('auto'),xlabel('x-position'),ylabel('y-position')
title('Zig-zag test')






