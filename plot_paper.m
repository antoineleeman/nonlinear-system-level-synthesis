
% File: polot_paper.m
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

load('data/05-Mar-2024_15_41_08__closed-loop.mat')

f=figure(1);

traj_nom.plot(1,"#0072BD");
traj_nom.plot_tubes(tubes + tubes_NL,"g",1,'Outer NL');
traj_nom.plot_tubes(tubes,'b',1,'LTV');
traj_noise.plot(1,'#D95319');

load('data/05-Mar-2024_15_49_05__nominal.mat');
traj_noise_nom.plot(1,'#B2BEB5');

l = traj_nom.add_constraints(m,1);
set(gcf,'units','centimeters','Position', [0 0 20 9]);
exportgraphics(gcf,'img/plot1.pdf','ContentType','vector');
saveas(gcf, 'fig1.png');

%%
clear all;
close all;
clc;

% Import CasADi toolbox
import casadi.*

load('data/05-Mar-2024_15_41_08__closed-loop.mat')

f = figure(3);
hold on;
s = subplot(1,2,1);
traj_nom.plot_tubes_phase(tubes+tubes_NL,'k','This paper (a)',traj_noise,s);
ylabel('\omega_2');
xlabel('\omega_3');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


clear all;
load('data/05-Mar-2024_15_41_50__open-loop.mat')

hold on;
s = subplot(1,2,2);
traj_nom.plot_tubes_phase(tubes+tubes_NL,'k','Open-loop (b)',traj_noise,s);
ylabel('\omega_2');
xlabel('\omega_3');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

ax = axes('Position', [0 0 1 1], 'Visible', 'off');
custom_colors = turbo(11);
colormap(turbo(11));
cb=colorbar;

mainAx = gca;
mainAx.Position = [0.3 0.2 0.6 0.65];
mainAx.CLim = [1 10];
mainAx.XAxis.Visible = 'off';
mainAx.YAxis.Visible = 'off';
ax = axes('Visible', 'off');
set(gcf,'units','centimeters','Position', [0 0 10 5]);

cb.Position = [0.9428    0.2007    0.0265    0.5479];
cb.Title.String = 'k';
cb.FontSize = 8;
exportgraphics(gcf,'img/plot2.pdf','ContentType','vector');