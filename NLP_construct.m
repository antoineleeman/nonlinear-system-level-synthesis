% File: NLP_construct.m
% Author: Antoine Leeman (aleeman(at)ethz(dot)ch)
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
classdef NLP_construct
    
    properties
        y; % vector of decision variables
        g; % vector of constraints
        f; % objection function
        x0; % initial guess
        n_y; % size of the vector y
        n_eq; % number of equality constraints within g
        n_ineq; % number of inequality constraints within g
        lbg; % vector of lower bounds for g
        ubg; % vector of upper bounds for g
        tol; % convergence tolerance
        inexact_jac; % (boolean) use of the ineaxact variant of the sequential quadratic programmig
        gamma;
    end
    
    methods
        function obj = NLP_construct(y,g,f,x0,n_y,n_eq, n_ineq)
            obj.y = y;
            obj.g = g;
            obj.f = f;
            obj.x0 = x0;
            obj.n_y = n_y;
            obj.lbg = [zeros(n_eq,1);-inf(n_ineq,1)];
            obj.ubg = zeros(n_eq +n_ineq,1 );
            obj.n_eq = n_eq;
            obj.n_ineq = n_ineq;
            obj.tol = 1e-6;
            obj.gamma = 1e-2;
        end
        
        function [res,output] = solve_nlp(obj)
            nlp = struct('x',obj.y, 'f',obj.f, 'g',obj.g);
            options.verbose =0;      

            options.ipopt.tol = obj.tol;
            options.ipopt.acceptable_tol = obj.tol;
            solver = casadi.nlpsol('solver', 'ipopt', nlp,options);

            res = solver('x0',obj.x0,'lbg',obj.lbg,'ubg',obj.ubg);
            output.t_proc_solver = solver.stats.t_proc_total;
            disp('converged');
            y = obj.y;
            g = obj.g;
            g_fun = casadi.Function('g_fun',{y},{g});
        end
        
        function [res_qp, output] = solve_sqp(obj,inexact_jac,n_eq_1)
            import casadi.*
            qp_solver = 'gurobi';
            gamma = obj.gamma;
            
            y = obj.y;
            g = obj.g;
            f = obj.f;
            
            y_bar = MX.sym('y_bar',length(y),1);
            Dy = MX.sym('Dy',length(y),1);
            
            dGdy = jacobian(g,y);          
            dGdy_exact = jacobian(g,y);
            if inexact_jac
                dGdy(n_eq_1(1):n_eq_1(2), 1:n_eq_1(3)) = casadi.DM.zeros(n_eq_1(2)-n_eq_1(1)+1,n_eq_1(3) ); % remove part of the Jacobians, according to Eq. (36)
            end
            dGdy_fun = casadi.Function('dGdy_fun',{y},{dGdy});
            dGdy_fun_exact = casadi.Function('dGdy_fun_exact',{y},{dGdy_exact});

            g_fun = casadi.Function('g_fun',{y},{g});           
            g_lin = g_fun(y_bar) + dGdy_fun(y_bar) * Dy;

            H_reg = gamma*casadi.DM.eye(length(y));
            H = evalf(hessian(f,y));
            H_tilde = H+H_reg;
            
            f_fun = casadi.Function('f_fun',{y},{f});
            dfdy = jacobian(f,y);
            dfdy_fun = casadi.Function('dfdy_fun',{y},{dfdy});
            
            f_qp = dfdy_fun(y_bar)* Dy + 0.5*Dy'*H_tilde*Dy;
            qp = struct('x',Dy, 'f',f_qp, 'g',g_lin,'p',y_bar);            

            
            switch qp_solver
                case 'ipopt'                   
                    options.ipopt.tol = 1e-10;
                    solver = casadi.nlpsol('solver', 'ipopt', qp, options);
                case 'gurobi'
                    options.gurobi.BarConvTol = 1e-9;
                    options.gurobi.FeasibilityTol = 1e-9;                    
                    solver = casadi.qpsol('solver', 'gurobi', qp, options);
            end
                                        
            iter_MAX = 50;
            y_current = zeros(length(y),iter_MAX+1);
            nu_current = zeros(length(g),iter_MAX+1);
            Dy = zeros(length(y),iter_MAX);
            p_feasibility = zeros(length(g),iter_MAX);  
            d_feasibility = zeros(length(g),iter_MAX);  
            
            primal_dual_step = zeros(1,iter_MAX);
            y_current(:,1) = zeros(length(y),1);

            converged = false;            
            for i = 1:iter_MAX

                res_qp = solver('lbg',obj.lbg,'ubg',obj.ubg,'p',y_current(:,i));
                
                y_current(:,i+1) = y_current(:,i)+ full(res_qp.x);
                Dy(:,i) =full(res_qp.x);
                                
                vec_g = g_fun(y_current(:,i));
                vec_g_eq = vec_g(1:obj.n_eq);
                vec_g_ineq = vec_g(obj.n_eq+1:end);
                vec_nu =  nu_current(:,i);
                vec_nu_ineq = vec_nu(obj.n_eq+1:end);

                p_feasibility(:,i) = full([vec_g_eq;vec_g_ineq]);
                
                primal_dual_step(i) = max([p_feasibility(:,i); d_feasibility(:,i)]);
                switch qp_solver
                    case 'gurobi'
                        t_proc_solver(i) =  solver.stats.t_wall_solver;
                    case 'ipopt'
                        t_proc_solver(i) = solver.stats.t_proc_total;
                end
                
                fprintf('Iteration %d, Primal-dual step: %f, Dual available %d\n', i, primal_dual_step(i), dual_available);
                if primal_dual_step(i) <= 1e-6
                    converged = true;
                    disp('Converged!')
                    res_qp.x = y_current(:,i);  
                    fprintf('Total solve time %f, Optimal cost: %f\n', sum(t_proc_solver), full(f_fun(res_qp.x)) );
                    break;
                end
            end
            
            output.p_feasibility = p_feasibility;
            output.t_proc_solver = t_proc_solver;
            output.converged = converged;
            output.primal_dual_step = primal_dual_step;
            output.optimal_cost = full(f_fun(res_qp.x));
                    
        end
        
    end
end

