% File: exp02_nominal_trajopt.m
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
N = m.N; % number of control intervals
nx = m.nx;
nu = m.nu;
ni = m.ni; % number of constraints
nw = m.nw; % size of the noise matrix E

% Create SLS object
sls = SLS(m);
sls.open_loop = false;

% Initialization of the solver parameters
solver_sqp = false; %false is ipopt
alpha = sls.alpha;

% Create CasADi symbolic variables
Z = MX.sym('state',nx,N+1);
V = MX.sym('input',nu,N+1);

% Initialization of the NLP
n_eq = 0; %number of equality constraints
n_ineq = 0; %number of inequality constraints

g_eq = []; %vector with all the equality constraints
g_ineq = []; %vector with all the inequality constraints


%% Initialization of the nominal trajectory optimization

% Get the nominal objective function (30a)
f = sls.getObjectiveNominal(m,Z,V);

% Dynamic constraint
[n_dyn,g_dyn] = sls.getConstraintsDynamic(m,Z,V);

n_eq = n_eq + n_dyn;
g_eq = [g_eq ; g_dyn];

% Nominal constraint satisfaction
[n_cons, g_cons, slack_cons, n_slack_cons] = sls.getConstraints(m,Z,V); 
n_ineq = n_ineq + n_cons;
g_ineq = [g_ineq ; g_cons];


% Form a vector with all the variables
[y_nom, n_y_nom] = sls.getVariablesNominal(Z,V);

y = [y_nom];
g = [g_eq;g_ineq]; % vector with all the constraints, with 0<=g_eq<=0 and -Inf <=g_eq<=0

f = f+ alpha*(y'*y); % modification of the objective function to ensure the Hessian is psd


%% SOLVE NLP
n_y = n_y_nom;
x0 = zeros(n_y,1); % initilisation of the SQP

nlp_cons = NLP_construct(y,g,f,x0,n_y,n_eq,n_ineq);

[res,output] = nlp_cons.solve_nlp(); %solve with IPOPT

z_nominal = res.x(1:nx*(N+1));
v_nominal = res.x(nx*(N+1)+1:end);

z = full(sls.vecToNominalState(z_nominal));
v = full(sls.vecToNominalInput(v_nominal));


%% SIMULATION
w = [    1     1     1     1     1     1     1     1     1    -1
     1     1     1     1     1     1    -1    -1    -1    -1
     1     1     1     1     1     1     1    -1    -1     1];

W = m.E*w;
x_cl_nom = zeros(nx,m.N+1);
x_cl_nom(:,1) = m.x0;

for k = 1 : m.N
    x_cl_nom(:,k+1)= m.ddyn_nl(x_cl_nom(:,k),v(:,k)) +W(:,k);
end

traj_noise_nom = Trajectory(m,x_cl_nom,v);
traj_nominal = Trajectory(m,z,v);

save(getUniqueName('nominal.mat'));
