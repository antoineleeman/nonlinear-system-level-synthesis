% File: main.m
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
addpath(genpath('../casadi-3.6.4_osx_arm64-matlab2018b')); % update with your casadi version!
%addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'))
addpath('util');


exp01_closed_loop_robust
exp02_nominal_trajopt
exp03_open_loop_robust

%todo: offline + multiplicative

plot_paper;