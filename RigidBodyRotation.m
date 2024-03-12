% File: RigidBodyRotation.m
% Author: Antoine Leeman (aleeman(at)ethz(dot)chh)
% Date: 05th March 2024
% License: MIT
% Reference:
%{
@article{leeman2023robust,
  title={Robust Nonlinear Optimal Control via System Level Synthesis},
  author={Leeman, Antoine P and K{\"o}hler, Johannes and Zanelli, Andrea and Bennani, Samir and Zeilinger, Melanie N},
  journal={arXiv preprint arXiv:2301.04943},
  year={2023}}
%}
% Link: https://arxiv.org/abs/2301.04943
% -----------------------------------------------------------------------------
%%
classdef RigidBodyRotation
    
    properties
        
        nx = 7; % number of state variables      
        nu = 3;% number of input variables
        ni = 6; % number of constraints        
        nw = 3;% size of the disturbance vector
        
        % initial and final states
        x0;
        xf;
        uf;
        % disturbance matrix
        E;
        
        % maximum torque and angular velocity values
        T_max;
        w_max;
        
        % input constraints matrices
        F_u;
        b_u;
        % state constraints matrices
        F_x;
        b_x;
        
        % number of time steps and total time
        N = 10;
        T = 10;
        dt;  % time step
        
        % cost matrices for state and input
        Q_cost;
        R_cost;
        
        mu; % constant which depends on the system (Eq. 9)
        I = 1*diag([5,2,1]); %inertia matrix
    end
    
    methods
            % constructor function
        function obj = RigidBodyRotation()
            
            q_0 = euler_to_quat([180 45 45]'); % convert initial orientation from Euler angles to quaternions
            w_0 = deg2rad([-1 -4.5 4.5]'); % convert initial angular velocity from degrees to radians
            q_f = euler_to_quat([0 0 0]'); % convert final orientation from Euler angles to quaternions
            w_f = deg2rad([0 0 0]'); % convert final angular velocity from degrees to radians
            
            % set initial state
            obj.x0 = [q_0;w_0];
            % set final state
            obj.xf = [q_f;w_f];
            obj.uf = zeros(3,1);
            
            % set maximum torque and angular velocity
            T_max= 0.1;
            w_max= 0.1;
            
            obj.T_max = T_max;
            obj.w_max = w_max;
            
            F_u = [eye(3); -eye(3)];
            b_u = T_max*[ones(6,1)];
            F_x = [[zeros(6,4)], [ eye(3);-eye(3)]];
            b_x = [w_max* ones(6,1)];
                        
            % set input constraints matrices
            obj.F_u = F_u;
            obj.b_u = b_u;
            % set state  constraints matrices
            obj.F_x = F_x;
            obj.b_x = b_x;

            obj.ni = length(b_u) + length(b_x); % total number of constraints

            % set cost matrices for state and input
            obj.Q_cost = 0.7*eye(obj.nx);
            obj.R_cost = eye(obj.nu);

            % set disturbance matrix
            E = zeros(obj.nw,obj.nx);
            E(1,5)= 0.005;E(2,6)= 0.005;E(3,7)= 0.005;
            obj.dt = obj.T/obj.N;
            obj.E = E'*obj.dt;

            obj.mu = [3.6993    3.7035    3.7170    3.6350    0.6492    4.6086    5.6358]; % To get those numerical values, we used m.compute_mu(10000)
        end
        
   
        function x_p = ddyn_nl(obj,x,u,integrator) %discretization of the dynamical system
            if nargin < 4
                integrator = 'rk4';
            end
            h = obj.dt;
            switch integrator
                case 'single'
                    x_p = x + h*ode_nl(obj,x,u);
                case 'multi'
                    step = 10;
                    for i = 1:step
                        x = x + h/step*ode_nl(obj,x,u);
                    end
                    x_p = x;
                case 'rk4'
                    k_1 = ode_nl(obj,x,u);
                    k_2 = ode_nl(obj,x+0.5*h*k_1,u);
                    k_3 = ode_nl(obj,x+0.5*h*k_2,u);
                    k_4 = ode_nl(obj,x+h*k_3,u);
                    x_p = x + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
                otherwise
                    error('unrecognised integrator');
            end
        end

        function dt = ode_nl(obj,x,u) % equation of motion of the dynamical system in continuous time (x_dot = ode(x,u) )
            I = obj.I;
            q = x(1:4);
            w = x(5:7);
            dt = [0.5 * quat_mut([0 w']',q);
                I\(u - cross(w,I*w))];
        end


        function A = A(obj,x,u) %state matric of the discrete time linearized dynamics
            import casadi.*
            x_fun = SX.sym('x',obj.nx);
            u_fun = SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            A = jacobian(obj.ddyn_nl(x_fun,u_fun), x_fun);
            A_fun = casadi.Function('A_fun',{var_fun},{A});
            A = A_fun([x;u]);
        end
        
        
        function B = B(obj,x,u) %input matrix of the discrete time linearized dynamics
            import casadi.*
            x_fun = SX.sym('x',obj.nx);
            u_fun = SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            B = jacobian(obj.ddyn_nl(x_fun,u_fun), u_fun);
            B_fun = casadi.Function('B_fun',{var_fun},{B});
            B = B_fun([x;u]);
        end

        function [max_mu] = compute_mu(obj,n_points)
            % estimation of the value of mu, via sampling
            rng(0,'twister');
            tic
            M = [ones(1,4), ones(1,3)*obj.w_max, ones(1,3)*obj.T_max];
            max_mu = zeros(1,obj.nx);
            parfor i = 1:n_points % assume symmetrical constraints
                eval = M.*(2*rand(1,10)-1);
                eval(1:4) = eval(1:4)/norm(eval(1:4)); % (!) non-uniform sampling
                max_mu = max(max_mu,obj.eval_mu(eval));
            end
            toc;
            
        end

        function mu = eval_mu(obj,xu)
            import casadi.*
            x_fun = SX.sym('x',obj.nx);
            u_fun = SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            H = jacobian(jacobian(obj.ddyn_nl(x_fun,u_fun,'rk4'), var_fun), var_fun);
            H_fun = casadi.Function('H_fun',{var_fun},{H});
            
            H = permute(reshape(full(H_fun(xu)),[obj.nx,obj.nx+obj.nu,obj.nx+obj.nu]),[3,2,1]);
            d = size(H);
            for i = 1:d(3)
                mu(i) = 0.5*sum(sum(abs(H(:,:,i))));
            end
        end
        
    end
end

