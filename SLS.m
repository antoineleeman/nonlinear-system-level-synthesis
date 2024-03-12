% File: SLS.m
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
classdef SLS
    properties
        N; % predictive horizon
        nx; % number of states
        nu; % number of input
        open_loop; % SLS_problem formulated in open_loop or closed_loop (boolean)
        linearized;% SLS_problem formulated for the original nonlinear dynamics or its linearization (boolean)
        alpha;
    end
    methods
        function obj = SLS(m)
            obj.N = m.N;
            obj.nx = m.nx;
            obj.nu = m.nu;
            obj.alpha = 1e-2;

        end
        
        function J = getObjectiveNominal(obj,m,Z,V)
            J =0;
            N = obj.N;
            for i=1:N
                J = J+V(:,i)'*m.R_cost*V(:,i) + (Z(:,i)-m.xf)'*m.Q_cost*(Z(:,i)-m.xf);
            end
            J = J+ (Z(:,N+1)-m.xf)'*m.Q_cost*(Z(:,N+1)-m.xf);
        end
        function [n,g] = getConstraintsDynamic(obj,m,Z,V)           
            g = [Z(:,1) - m.x0];
            n = obj.nx;
            
            for i=1:obj.N
                g = [g; Z(:,i+1) - m.ddyn_nl(Z(:,i),V(:,i))];
            end
            n = n+ obj.N*obj.nx;
        end

        function [n_ineq,g_ineq,var_cons,n_var_cons] = getRobustConstraints(obj,m,Z,V, v3_Phi_x,v3_Phi_u, delta)
            N = obj.N;
            C_x = m.F_x;
            D_x = m.b_x;
            
            C_u = m.F_u;
            D_u = m.b_u;
            
            var_cons = [];
            n_var_cons =0;
            n_ineq = 0;
            g_ineq = [];
            for k=1:N
                Phi_k = SLS.Phi_line_k(k, v3_Phi_x);
                CPhi_x = C_x * Phi_k ;
                E_kron = kron(eye(k),[zeros(7,m.nx-m.nw), m.E]);

                delta_k = delta(1:k).^2;
                delta_kron = kron(diag(delta_k), diag(m.mu)) + E_kron;
                [one_norm_2, n,g, var, n_var_slack] = SLS.vec_one_norm(CPhi_x* delta_kron);
                g_ineq = [g_ineq;g];
                
                n_ineq = n_ineq+n;
                var_cons = [var_cons;var;one_norm_2];
                n_var_cons = n_var_cons + n_var_slack;
                g_ineq = [g_ineq; C_x*Z(:,k+1)+one_norm_2- D_x];
                
                n_ineq = n_ineq + length(D_x);
            end
            for k=1:N
                Phi_k = SLS.Phi_line_k(k, v3_Phi_u);
                CPhi_u = C_u * Phi_k;
                E_kron = kron(eye(k),[zeros(7,m.nx-m.nw), m.E]);
                
                delta_k = delta(1:k).^2;
                delta_kron = kron(diag(delta_k), diag(m.mu)) + E_kron;
                [one_norm_2, n,g, var, n_var_slack] = SLS.vec_one_norm(CPhi_u* delta_kron);
                g_ineq = [g_ineq;g];
                n_ineq = n_ineq+n;
                var_cons = [var_cons;var;one_norm_2];
                n_var_cons = n_var_cons + n_var_slack;
                
                g_ineq = [g_ineq; C_u*V(:,k+1)+ one_norm_2  - D_u];
                
                n_ineq = n_ineq + length(D_u);
            end
            g_ineq = [g_ineq; C_u*V(:,1) - D_u;C_x*Z(:,1) - D_x; ];
            
            n_ineq = n_ineq + m.ni;
            
        end

        function [n_ineq,g_ineq,var_cons,n_var_cons] = getConstraints(obj,m,Z,V)
            N = obj.N;
            C_x = m.F_x;
            D_x = m.b_x;

            C_u = m.F_u;
            D_u = m.b_u;

            var_cons = [];
            n_var_cons =0;
            n_ineq = 0;
            g_ineq = [];
            for k=1:N
                g_ineq = [g_ineq; C_x*Z(:,k+1)- D_x];
                n_ineq = n_ineq + length(D_x);
            end
            for k=1:N

                g_ineq = [g_ineq; C_u*V(:,k+1)- D_u];

                n_ineq = n_ineq + length(D_u);
            end
            g_ineq = [g_ineq; C_u*V(:,1) - D_u;C_x*Z(:,1) - D_x; ];

            n_ineq = n_ineq + m.ni;

        end


        function [n,g] = getLTVMapConstraints(obj, m,Z,V, Phi_x, Phi_u)
            N = obj.N;
            nx = obj.nx;
            nu = obj.nu;
            
            R = obj.v3_to_R(Phi_x);
            M = obj.v3_to_M(Phi_u);
            I = kron(eye(nx),eye(N));
            
            AB = obj.buildLTVMap(m,Z,V);
            
            g_full = [reshape(AB*[R;M] - I, [N^2*nx^2,1])];
            [n,g] = obj.reducedMap(g_full);
        end
        function [n_ineq,g_ineq,n_var_norm,var_norm] = getConstraintTube(obj,m,delta, v3_Phi_x, v3_Phi_u)
            g_ineq = [];
            n_ineq = 0;
            n_var_norm = 0;
            var_norm = [];
            
            for k = 1:obj.N-1
                Phix_k = SLS.Phi_line_k(k, v3_Phi_x);
                Phiu_k = SLS.Phi_line_k(k, v3_Phi_u);
                Phi_k = [Phix_k;Phiu_k];
                
                E_kron = kron(eye(k),m.E);
                [ inf_norm, n_ineq_norm,g_ineq_norm, var_slack_norm,n_var_slack_norm ] = obj.mat_inf_norm(Phi_k* E_kron);
                n_var_norm = n_var_norm + n_var_slack_norm;
                var_norm = [var_norm;var_slack_norm;inf_norm];
                n_ineq = n_ineq + n_ineq_norm;
                g_ineq = [g_ineq; g_ineq_norm];
                
                delta_k = delta(1:k).^2;
                delta_kron = kron(diag(delta_k), diag(m.mu));
                [ inf_norm_2, n_ineq_norm,g_ineq_norm, var_slack_norm,n_var_slack_norm ] = obj.mat_inf_norm(Phi_k* delta_kron);
                n_var_norm = n_var_norm + n_var_slack_norm;
                var_norm = [var_norm;var_slack_norm;inf_norm_2];
                n_ineq = n_ineq + n_ineq_norm;
                g_ineq = [g_ineq; g_ineq_norm];
                
                g_ineq = [g_ineq; inf_norm + inf_norm_2- delta(k+1)];
                n_ineq = n_ineq+1;
            end
        end
                
        function AB = buildLTVMap(obj,m,Z,V)
            % Construction of the matrix [I-ZA, -ZB] line by line.
            nx = obj.nx;
            nu = obj.nu;
            N = obj.N;
            
            AB = sparse(nx, N*(nx+nu));
            AB(1:nx,1:nx) = eye(nx); % first line
            
            for i =1:N-1 % loop between the second line until the last.
                A = m.A(Z(:,i+1),V(:,i+1));
                B = m.B(Z(:,i+1),V(:,i+1));
                mat_1 = sparse(nx,(i-1)*nx);
                mat_2 = [-A, eye(nx)];
                mat_3 =  sparse(nx,(N-2)*nx + (i-1)*(nu-nx));
                mat_4 = -B;
                mat_5 = sparse(nx, (N-i)*nu );
                AB_ = [mat_1, mat_2, mat_3,mat_4,mat_5];
                AB = [AB;AB_];
            end
        end       
        function [n,g] = reducedMap(obj,g)
            % Remove all the unnecessary equality constraints (i.e., 0=0)
            % from the System Level Parametrization
            N = obj.N;
            nx = obj.nx;
            v_NZ = find(reshape(kron(tril(ones(N),0),ones(nx)), [N^2*nx^2,1])); % remove all trivial equality constraints
            g = g(v_NZ);
            n = N*(N+1)/2*nx^2;
        end
        
        function [z_sol,v_sol,Phi_x_sol,Phi_u_sol,tube_sol,f_sol]= solConvertVecMat(obj,m,y_sol, n_y_nom, n_y_contr, n_y_tube)
            N = obj.N;
            nx = obj.nx;
            n=0;
            z_sol_v = full(y_sol(n+1:(N+1)*nx));
            v_sol_v = full(y_sol((N+1)*nx+1:n_y_nom));
            n = n+n_y_nom;
            Phi_sol_v = full(y_sol(n+1 : n +n_y_contr));
            n = n+n_y_contr;
            tube_sol = full(y_sol(n+1 : n+n_y_tube));
            
            z_sol = obj.vecToNominalState(z_sol_v);
            v_sol = obj.vecToNominalInput(v_sol_v);            
            f_sol = obj.getObjectiveNominal(m,z_sol,v_sol);
            
            [Phi_x_sol,Phi_u_sol] = obj.vecToResponse(Phi_sol_v);
            
        end        
        function [y,n] = getVariablesNominal(obj,Z,V)
            N = obj.N;
            nx = obj.nx;
            nu = obj.nu;
            y = [reshape(Z,[(N+1)*nx,1]); reshape(V,[(N+1)*nu,1])];
            n = (N+1)*(nx+nu);
        end       
        function [y,n] = getVariablesTube(obj,delta)
            N = obj.N;
            y = reshape(delta,[N,1]);
            n = N;
        end
        function [y,n] = getVariablesResponses(obj,Phi_x,Phi_u)
            y = [];
            N = obj.N;
            nu = obj.nu;
            nx = obj.nx;
            n=0;
            if ~obj.open_loop
                for i = 1:(N+1)*N/2
                    y = [y; reshape(Phi_u{i},[nu*nx,1])];
                end
                n = n + (N+1)*N/2 * nu*nx;
            end
            for i = 1:(N+1)*N/2
                y = [y;reshape(Phi_x{i},[nx^2,1])];
            end
            n = n + (N+1)*N/2*nx^2;
        end
        
        function Z = vecToNominalState(obj,y)%
            nx = obj.nx;
            N = obj.N;
            Z = reshape(y,[nx,N+1]);
        end        
        function V = vecToNominalInput(obj,y)%
            nu = obj.nu;
            N = obj.N;
            V = reshape(y,[nu,N+1]);
        end
        function [Phi_x,Phi_u] = vecToResponse(obj,y)%
            nu = obj.nu;
            N = obj.N;
            nx = obj.nx;
            k = 1;
            if obj.open_loop
                Phi_u = cellmat(N*(N+1)/2,1,nu,nx);
            else
                for n =1: (N+1)*N/2
                    Phi_u{n} = reshape(y(k:k+nu*nx-1),[nu,nx]);
                    k = k + nu*nx;
                end
            end
            for n =1: (N+1)*N/2
                Phi_x{n} = reshape(y(k:k+nx*nx-1),[nx,nx]);
                k = k + nx^2;
            end
        end
 
        function R = v3_to_R(obj,v3_Phi_x)%
            R = [];
            n=1;
            for k = 1:obj.N
                R_= [];
                for j=1:k
                    R_ = [R_, v3_Phi_x{n}];
                    n = n+1;
                end
                for j=k+1:obj.N
                    R_ = [R_, sparse(obj.nx,obj.nx)];
                end
                R = [R;R_];
            end
            
        end
        function M = v3_to_M(obj,v3_Phi_u)%
            M = [];
            n=1;
            for k = 1:obj.N
                M_ = [];
                for j=1:k
                    M_ = [M_, v3_Phi_u{n}];
                    n = n+1;
                end
                for j=k+1:obj.N
                    M_ = [M_, sparse(obj.nu,obj.nx)];
                end
                M = [M;M_];
            end
        end


    end
        methods(Static)
        function L = Phi_line_k(k, v3_Phi)%
            %k is the number of the line, starting from 1
            L = [];
            last = k*(k+1)/2;
            first = k*(k-1)/2+1;
            for j=first:last
                L = [L, v3_Phi{j}];
            end
        end

        function [one_norm_slack, n,g_ineq, slack, n_var] = one_norm(v)%
            
            import casadi.*
            n_slack = length(v);
            slack = MX.sym('slack_vec',n_slack); % 1-norm always computed on vertical vectors
            one_norm_slack = MX.sym('slack',1);
            %g_ineq = [ v - slack; -v - slack; sum(slack) - one_norm_slack;epsilon-one_norm_slack];
            g_ineq = [ v - slack; -v - slack; sum(slack) - one_norm_slack];
            
            n = 2*n_slack +1;%+1;%%
            n_var = n_slack +1;
        end        
        function [one_norm_line, n_ineq,g_ineq, var_slack,n_var_slack ] = vec_one_norm(M)%
            var_slack = [];
            n_var_slack = 0;
            n_ineq = 0;
            g_ineq = [];
            dim = size(M);
            one_norm_line = [];
            for i = 1 : dim(1)
                [one_norm,n,g, var, n_var] = SLS.one_norm(M(i,:)');
                one_norm_line = [one_norm_line; one_norm];
                n_ineq = n_ineq+n;
                g_ineq = [g_ineq;g];
                var_slack = [var_slack;var];
                n_var_slack  = n_var_slack + n_var;
            end
        end        
        function [inf_norm, n_ineq,g_ineq, var_slack,n_var_slack ] = mat_inf_norm(M)%
            import casadi.*
            inf_norm = MX.sym('slack_tube',1);
            
            var_slack = [];
            n_var_slack = 1;
            n_ineq = 0;
            g_ineq = [];
            dim = size(M);
            for i = 1 : dim(1)
                [one_norm,n,g, var, n_var] = SLS.one_norm(M(i,:)');
                n_ineq = n_ineq+n;
                g_ineq = [g_ineq;g];
                var_slack = [var_slack;var;one_norm];
                n_var_slack  = n_var_slack + n_var;
                
                g_ineq = [g_ineq;one_norm - inf_norm];
                n_ineq = n_ineq + 1;
            end
        end
   
    end
end

