classdef Heli < handle
    % This class provides an object of the Bi-Rotor Planar Helicopter. 
    % Various controllers are provided to control the system:
    %
    % 1. Linear Quadratic Regulator
    %
    % 2. Linear Model Reference Adaptive Controller
    %
    % 3. Non-linear Model Reference Adaptive Controller
    %
    
    properties
        % All units are SI
        m; % mass of the rotor
        d; % mass flow rate
        g; % gravity 
        r; % distance from center of gravity
        J; % moment of inertia
    end
    
    methods
        function self = Heli(m,d,g,r,J)
            if nargin > 0
                self.m = m;
                self.d = d;
                self.g = g;
                self.r = r;
                self.J = J;
            else
                % These values will be used if no value is provided to the
                % class
                self.m = 6; % kg
                self.d = 0.1; % kg/sec
                self.g = 9.8; % m/sec^2
                self.r = 0.25; % m
                self.J = 0.1425; % kg m^2
            end
        end
        
        function out = linear_controller(self, t, x,A,B,K)
            % t = current time
            % x = current state
            % Q = state weight matrix
            % R = input weight matrix

            u = -K*x;
            % ode at current time step
            out = A*x + B*u;
        end
        
        function out = mrac_controller(self,t, x,A,B,Am,Bm,Q,R)
            kx = reshape(x(13:24),6,2);
            kr = reshape(x(25:28),2,2);
            ref = [3;0];
           
            u = kx'*x(1:6) + kr'*ref;

            gamma_x = 1*eye(6);
            gamma_r = 1*eye(2);
            e = x(1:6) - x(7:12);
            [P,~,~] = care(A, B, Q, R);
            kxdot = -gamma_x*x(1:6)*e'*P*B;
            krdot = -gamma_r*ref*e'*P*B;

            %Helicopter Model
            dXdt = A*x(1:6) + B*u;
            dXmdt = Am*x(7:12) + Bm*ref;
            out = [dXdt; dXmdt; reshape(kxdot,numel(kxdot),1);reshape(krdot,numel(krdot),1)];
        end
        
        function out = nl_mrac_controller(self, t, x, A, B, Am, Bm, Q, R, beta_nl)
            kx = reshape(x(13:48),6,6);
            kr = reshape(x(49:84),6,6);
            beta = reshape(x(85:120),6,6);
            ref = [1.5;1.5;0;0;0;0];
            u = kx'*x(1:6) + kr'*ref;
            phi = [ 1; x(6); 0; x(6)^2; 0; 0];
            
            gamma_x = 50*eye(6);
            gamma_r = 50*eye(6);
            gamma_b = 10*eye(6);
            
            e = x(1:6) - x(7:12);
            P = care(A, B, Q, R);
            kxdot = -gamma_x*x(1:6)*e'*P*B;
            krdot = -gamma_r*ref*e'*P*B;
            b_dot = -gamma_b*phi*e'*P*B;
            %Helicopter Model
            dXdt = A*x(1:6) + B*(u + beta_nl'*phi - beta'*phi);
            dXmdt = Am*x(7:12) + Bm*ref;
            out = [dXdt; dXmdt; reshape(kxdot,numel(kxdot),1);reshape(krdot,numel(krdot),1);reshape(b_dot,numel(b_dot),1)];
        end
    end
end

