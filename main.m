% Playground for the project

% Object of the helicopter controller
helicopter = Heli(m,d,g,r,J);

% System specifications
d = 0.1; m = 6; r = 0.25; J = 0.1425; g = 9.8; f = [0.5*g;0.5*g];

% Initial conditions and time
tspan = [0, 30];
x0 = [5 1 5 1 pi/4 pi/2];

% Weighing matrices for quadratic solver
Q = eye(6);
R = eye(2);
R_nl = eye(6);
% Thrusts
fr = [norm([f(1)*sin(x0(5)), f(1)*cos(x0(5)) - g]), norm([f(2)*sin(-x0(5)), f(2)*cos(x0(5)) - g])];
            
% System State matrices
A = [0 1 0 0 0 0; 0 -d/m 0 0 -(fr(1)+fr(2))/m 0; 0 0 0 1 0 0; 0 0 0 -d/m 0 0; 0 0 0 0 0 1; 0 0 0 0 0 0];
B = [0 0; 0 0; 0 0; 1/m 1/m; 0 0; r/J -r/J];

A_nl = [0 1 0 0 0 0;0 -d/m 0 0 0 0;0 0 0 1 0 0;0 0 0 -d/m 0 0;0 0 0 0 0 1;0 0 0 0 0 0];
B_nl = [0 0 0 0 0 0;0 cos(x(5)) 0 -sin(x(5)) 0 0;0 0 0 0 0 0;-1 sin(x(5)) 0 cos(x(5)) 0 0;0 0 0 0 0 0;0 0 0 0 0 1];
beta_nl = [ g   0       0   0   0   0;
            0   -d/m    0   0   0   0;
            0   0      -1   0   0   0;
            0   0       0   0   0   0;
            0   0       0   0   0   0;
            0   0       0   0   0   0;]';

% Solving Ricatti to find optimal gain
[~,~,K] = care(A,B,Q,R);

[~,~,K_nl] = care(A_nl,B_nl,Q,R_nl);
% K_s = [-3.8644 1.0268 -3.6296 -26.0352 2.7563 -5.0786; -1.0311 4.3841 -1.5035 -25.9090 -2.7249 7.5734];

% Model State matrices
Am = A-B*K;
Bm = B;

% Non-linear Model State matrice
Am_nl = A_nl-B_nl*K_nl;
Bm_nl = B_nl;

% kr_s = [-0.0608 0.0000; 0.4546 1.0000];
kr_star = (pinv(B)*Bm);  
kr_star_nl = (pinv(B_nl)*Bm_nl);  
% Intital conditions for adaptive controller
beta0 = ones([6,6]);
x0_adapt = [x0'; 4; 0; 4; 0; pi/6; pi/3; reshape(K',numel(K'),1); reshape(kr_star',numel(kr_star'),1)];
x0_adapt_nl = [x0'; 4; 0; 4; 0; pi/6; pi/3; reshape(K_nl',numel(K_nl'),1); reshape(kr_star_nl',numel(kr_star_nl'),1); reshape(beta0',numel(beta0'),1)];

% Differential system solver
% [t,x] = ode45(@(t,x) helicopter.mrac_controller(t, x, A, B, Am, Bm, Q, R), tspan, x0_adapt);

% [t,x] = ode45(@(t,x) helicopter.linear_controller(t, x, A, B, K), tspan, x0);

[t,x] = ode45(@(t,x) helicopter.nl_mrac_controller(t, x, A_nl, B_nl, Am_nl, Bm_nl, Q, R_nl, beta_nl), tspan, x0_adapt_nl);

%%%%%%%%%%%%%%% TESTING POST TRANSIENT GAINS %%%%%%%%%%%%%
% Since tspan is 50, which is very high; we assume that at t = t_n/2 we will
% have non-transient behaviour. Hence we take the gains from that time.

out_mid = x(floor(end/2), :);
K_pt = reshape(x(13:48),6,6);
% 
kr_pt = reshape(x(49:84),6,6);
%
beta_pt = reshape(out_mid(85:120),6,6);

x0_adapt_nl = [x0'; 4; 0; 4; 0; pi/6; pi/3; reshape(K_pt',numel(K_pt'),1); reshape(kr_pt',numel(kr_pt'),1);reshape(beta_pt',numel(beta_pt'),1)];

% Differential system solver
[t_pt,x_pt] = ode45(@(t,x) helicopter.nl_mrac_controller(t, x, A_nl, B_nl, Am_nl, Bm_nl, Q, R_nl,beta_nl), tspan, x0_adapt_nl);

% fr = zeros(length(t),2);
% 
% for temp = 1:length(t)
%     theta = wrapToPi(x(temp,5));
%     fr(temp,:) = [norm([f(1)*sin(theta), f(1)*cos(theta) - g]), norm([f(2)*sin(-theta), f(2)*cos(theta) - g])];
% end

% figure(1);
%     subplot(3,1,1);
%     hold on;
%     plot(t,x(:,1));
%     xlabel('t');
%     ylabel('x(t)');
%     hold off;
%     grid on;
%     
%     subplot(3,1,2);
%     hold on;
%     plot(t,x(:,3));
%     xlabel('t');
%     ylabel('y(t)');
%     hold off;
%     grid on;
%     
%     subplot(3,1,3);
%     hold on;
%     plot(t,x(:,5));
%     xlabel('t');
%     ylabel('theta(t)');
%     hold off;
%     grid on;
%     
% figure(2)
%   hold on;
%   plot(t,fr(:,1));
%   plot(t,fr(:,2));
%   title('Force applied');
%   legend('f_1','f_2');
%   hold off;
%   grid on;
% 
% figure(3);
%   subplot(1,1,1);
%   hold on
%   plot(x(:,1), x(:,3));
%   for temp = 1:floor(length(t)/10)
%      arrow([x(10*temp,1),x(10*temp,3)],[x(10*temp,1)+0.5*cos(x(10*temp,5)),x(10*temp,3)+0.5*sin(x(10*temp,5))]);
%   end
%   xlabel('x');
%   ylabel('y');
%   title("Trajectory and Orientation")
%   hold off;
%   grid on;

figure(1);
  subplot(3,1,1);
  hold on
  plot(t, x(:,1));
  plot(t_pt, x_pt(:,1));
  plot(t, x(:,7));
  xlabel('t');
  ylabel('x');
  legend('x','x_{post\_transient}','x_m');
  hold off
  grid on;
  
  subplot(3,1,2);
  hold on
  plot(t, x(:,3));
  plot(t_pt, x_pt(:,3));
  plot(t, x(:,9));
  xlabel('t');
  ylabel('y');
  legend('y','y_{post\_transient}','y_m');
  hold off
  grid on;
  
  subplot(3,1,3);
  hold on
  plot(t, x(:,5));
  plot(t_pt, x_pt(:,5));
  plot(t, x(:,11));
  xlabel('t');
  ylabel('{\theta}');
  legend('{\theta}','{\theta_{post\_transient}}','{\theta}_m');
  hold off
  grid on;
  
figure(2);
  subplot(1,1,1);
  hold on
  plot(t, x(:,1)-x(:,7));
  plot(t, x(:,2)-x(:,8));
  plot(t, x(:,3)-x(:,9));
  plot(t, x(:,4)-x(:,10));
  plot(t, x(:,5)-x(:,11));
  plot(t, x(:,6)-x(:,12));
  xlabel('t');
  ylabel('error');
  legend('e1','e2','e3','e4','e5','e6');
  hold off
  grid on;
  
figure(3);
  subplot(1,1,1);
  hold on
  plot(x(:,1), x(:,3));
  plot(x(:,7), x(:,9));
  plot(x_pt(:,1), x_pt(:,3));
  
%   for temp = 1:floor(length(t)/100)
%      arrow([x(100*temp,1),x(100*temp,3)],[x(100*temp,1)+0.5*cos(x(100*temp,5)),x(100*temp,3)+0.5*sin(x(100*temp,5))]);
%   end
%   
%   for temp = 1:floor(length(t)/500)
%      arrow([x(100*temp,7),x(100*temp,9)],[x(100*temp,7)+0.5*cos(x(100*temp,11)),x(100*temp,9)+0.5*sin(x(100*temp,11))]);
%   end
%   
%    for temp = 1:floor(length(t)/500)
%      arrow([x_pt(100*temp,7),x_pt(100*temp,9)],[x_pt(100*temp,7)+0.5*cos(x_pt(100*temp,11)),x_pt(100*temp,9)+0.5*sin(x_pt(100*temp,11))]);
%    end
  
  xlabel('x');
  ylabel('y');
  legend('system','model','system with post-transient gains');
  hold off
  grid on;
  
figure(4);
  subplot(2,1,1);
  hold on;
  
  for temp = 1:36
      plot(t,x(:,12+temp));
      xlabel('t');
      ylabel('kx gains');
  end
  
  hold off;
  
  subplot(2,1,2);
  hold on;
  for temp = 1:36
      plot(t,x(:,48+temp));
      xlabel('t');
      ylabel('kr gains');
  end
  hold off;