% File: exp03_open_loop_robust.m
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
clear all;
close all;
clc;

% Import CasADi toolbox
import casadi.*

% Create dynamics system Object
m = RigidBodyRotation();

m.N = 6;
m.T = 6;
N = m.N; % number of control intervals
nx = m.nx;
nu = m.nu;
ni = m.ni; % number of constraints
nw = m.nw; % size of the noise matrix E

% Create SLS object
sls = SLS(m);
sls.open_loop = true;

% Initialization of the solver parameters
inexact_jac=false;
use_jit = false;
solver_sqp = false; %false is ipopt

alpha = sls.alpha;

% Create CasADi symbolic variables
Z = MX.sym('state',nx,N+1);
V = MX.sym('input',nu,N+1);
Phi_x = MX.sym('state_response',nx,nx,N*(N+1)/2);
Phi_u = cellmat(N*(N+1)/2,1,nu,nx); %%OPEN-LOOP ROBUST!

delta = MX.sym('lin_tube',N);

% Initialization of the NLP
n_eq = 0; %number of equality constraints
n_ineq = 0; %number of inequality constraints

g_eq = []; %vector with all the equality constraints
g_ineq = []; %vector with all the inequality constraints


%% Initialization of the NLP (30)

% Get the nominal objective function (30a)
f = sls.getObjectiveNominal(m,Z,V);

% Dynamic constraint
[n_dyn,g_dyn] = sls.getConstraintsDynamic(m,Z,V); %Eq.(30c)

n_eq = n_eq + n_dyn;
g_eq = [g_eq ; g_dyn];

% Robust constraint satisfaction Eq. (30d)
[n_cons, g_cons, slack_cons, n_slack_cons] = sls.getRobustConstraints(m,Z,V, Phi_x,Phi_u, delta); 
n_ineq = n_ineq + n_cons;
g_ineq = [g_ineq ; g_cons];

% System Level Parametrization Eq.(30b)
[n_map, g_map] = sls.getLTVMapConstraints(m, Z,V, Phi_x,Phi_u);
n_eq = n_eq + n_map;
g_eq = [g_eq ; g_map];

[n_ineq_tube,g_ineq_tube, n_slack_tube, slack_tube] = sls.getConstraintTube(m,delta, Phi_x, Phi_u); %Eq. (30e)
n_ineq = n_ineq + n_ineq_tube;
g_ineq = [g_ineq ; g_ineq_tube];

% Form a vector with all the variables
[y_nom, n_y_nom] = sls.getVariablesNominal(Z,V);
[y_contr,n_y_contr] = sls.getVariablesResponses(Phi_x,Phi_u);

[y_tube, n_y_tube] = sls.getVariablesTube(delta);
var_slack = [slack_cons;slack_tube];

n_eq_1 = [n_eq-n_map+1,n_eq,n_y_nom]; % equality whose jacobians are neglected in the (inexact) SQP variant

y = [y_nom; y_contr;y_tube;var_slack];
g = [g_eq;g_ineq]; % vector with all the constraints, with 0<=g_eq<=0 and -Inf <=g_eq<=0

f = f+ alpha*(y'*y); % modification of the objective function to ensure the Hessian is psd

n_y = n_slack_cons+n_slack_tube+ n_y_tube+n_y_contr + n_y_nom;
x0 = zeros(n_y,1); % initilisation of the SQP

nlp_cons = NLP_construct(y,g,f,x0,n_y,n_eq,n_ineq);
nlp_cons.use_jit = use_jit;
nlp_cons.inexact_jac = inexact_jac;

%% Solution of the NLP
if ~solver_sqp
    [res,output] = nlp_cons.solve_nlp(); %solve with IPOPT
else
    if nlp_cons.inexact_jac
        [res, output] = nlp_cons.solve_sqp(lintrue,n_eq_1); % solve with inexact SQP variant
    else
        [res, output] = nlp_cons.solve_sqp(false); % solve with classic SQP
    end
end

%% Post-processing of the solution
[z_sol,v_sol,Phi_x_sol,Phi_u_sol,tube_sol,f_sol]= sls.solConvertVecMat(m,res.x, n_y_nom, n_y_contr, n_y_tube);

R_sol = sls.v3_to_R(Phi_x_sol);
M_sol = sls.v3_to_M(Phi_u_sol);

min_delta(1)=0;
for k = 1:sls.N-1
    Phix_k = sls.Phi_line_k(k, Phi_x_sol);
    Phiu_k = sls.Phi_line_k(k, Phi_u_sol);
    Phi_k = [Phix_k;Phiu_k];
    
    E_kron = kron(eye(k),[m.E]);

    inf_norm_1 = norm(Phi_k* E_kron,Inf);
    delta_k = min_delta(1:k).^2;
    delta_kron = kron(diag(delta_k), diag(m.mu));
    inf_norm_2 = norm(Phi_k* delta_kron,Inf);
    
    min_delta(k+1) = inf_norm_1+inf_norm_2;
end

tubes_x = reshape(vecnorm( kron(eye(m.N),[eye(m.nx);-eye(m.nx)]) * R_sol * kron(eye(m.N),m.E),1,2), [2*m.nx,m.N]);
tubes_x = [zeros(2*m.nx,1), tubes_x];

tubes_u = reshape(vecnorm( kron(eye(m.N),m.F_u) * M_sol * kron(eye(m.N),m.E),1,2), [2*m.nu,m.N]);
tubes_u = [zeros(2*m.nu,1), tubes_u];

delta_kron = kron(diag(min_delta.^2), diag(m.mu));
tubes_x_NL = reshape(vecnorm( kron(eye(m.N),[eye(m.nx);-eye(m.nx)]) * R_sol * delta_kron ,1,2 ), [2*m.nx,m.N]);
tubes_x_NL = [zeros(2*m.nx,1), tubes_x_NL];

tubes_u_NL = reshape(vecnorm( kron(eye(m.N),m.F_u) * M_sol * delta_kron ,1,2), [2*m.nu,m.N]);
tubes_u_NL = [zeros(2*m.nu,1), tubes_u_NL];


tubes = [tubes_x;tubes_u];

tubes_NL = [tubes_x_NL; tubes_u_NL];

K = M_sol/R_sol;

%% SIMULATION

C_x = m.F_x;
C_x = C_x(5,:);
Phi_k = sls.Phi_line_k(N, Phi_x_sol);
w = reshape(sign(C_x*Phi_k*kron(eye(N),m.E)),[m.nw,sls.N]);

W = [m.E]*w;

x_cl = zeros(nx,m.N+1);
x_cl(:,1) = m.x0;
for k = 1 : m.N
    Delta_x = reshape(x_cl(:,2:end) -z_sol(:,2:end), [N*nx,1]);
    u_cl = reshape([zeros(nu,1);K*Delta_x],[nu,N+1]) + v_sol;
    x_cl(:,k+1)= m.ddyn_nl(x_cl(:,k),u_cl(:,k)) +W(:,k);
end
traj_noise = Trajectory(m,x_cl,u_cl);
traj_nom = Trajectory(m,z_sol,v_sol);


% SAVE
save(getUniqueName('open-loop.mat'));
