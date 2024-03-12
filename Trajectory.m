% File: Trajectory.m
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
classdef Trajectory
    % Represents the states and inputs of a system over a certain time
    properties
        X; % states
        U; % input
        T; % time vector
        m;  % model object
        dt; % time step
        K;  % number of steps
    end
    
    methods
        function obj = Trajectory(m,X,U) % Constructor method
            obj.X = X;
            obj.U = U;
            obj.dt = m.dt;
            obj.T = 0:m.dt:m.dt*(m.N);
            obj.m = m;
            obj.K = length(obj.T);
        end
        
        function obj = plot(obj,fig_num,color) % use stackedplot !!
            
            figure(fig_num);
            subplot(3,3,[1,4]);
            hold on
            grid on
            plot(obj.T,obj.X(1:4,:)','.-','Color',color);
            xlabel('Step $k$','Interpreter','latex')
            %ylabel('Quaternion ${q}$','Interpreter','latex')
            set(gca,'fontsize',7);
            
            subplot(3,3,2);
            subtitle('Angular speed ${\omega}$','Interpreter','latex')
            hold on
            grid on

            %stackedplot(obj.T, obj.X(5:7,:)' )
            plot(obj.T,obj.X(5,:),'.-','Color',color);
            ylabel('${\omega}_1$','Interpreter','latex')

            set(gca,'fontsize',7);
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])

            subplot(3,3,5);
            hold on
            grid on
            plot(obj.T,obj.X(6,:),'.-','Color',color);
            ylabel('${\omega}_2$','Interpreter','latex')
            set(gca,'fontsize',7);
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])

            subplot(3,3,8);
            hold on
            grid on
            plot(obj.T,obj.X(7,:),'.-','Color',color);
            xlabel('Step $k$','Interpreter','latex')
            ylabel('${\omega}_3$','Interpreter','latex')
            
            set(gca,'fontsize',7);
            
            subplot(3,3,3);
            hold on
            grid on
            stairs(obj.T',[obj.U(1,1:end-1), obj.U(1,end-1)]','.-','Color',color)
            %xlabel('Step $k$','Interpreter','latex')
            subtitle('Control input ${u}$','Interpreter','latex');
            ylabel('${u}_1$','Interpreter','latex')
                        set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'fontsize',7);
            
            subplot(3,3,6);
            hold on
            grid on
            stairs(obj.T',[obj.U(2,1:end-1), obj.U(2,end-1)]','.-','Color',color)
            %xlabel('Step $k$','Interpreter','latex')
            set(gca,'fontsize',7);
            ylabel('${u}_2$','Interpreter','latex')
            
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            
            subplot(3,3,9);
            hold on
            grid on
            stairs(obj.T',[obj.U(3,1:end-1), obj.U(3,end-1)]','.-','Color',color)
            xlabel('Step $k$','Interpreter','latex')
            set(gca,'fontsize',7);
            ylabel('${u}_3$','Interpreter','latex')
            

        end
        
        
        function obj = plot_stacked(obj,fig_num,color)

            figure(fig_num);
            subplot(3,3,[1,4]);
            hold on
            grid on
            plot(obj.T,obj.X(1:4,:)','.-','Color',color);
            xlabel('Step $k$','Interpreter','latex')
            %ylabel('Quaternion ${q}$','Interpreter','latex')
            set(gca,'fontsize',7);
            
            subplot(3,3,[2,5,8]);
            subtitle('Angular speed ${\omega}$','Interpreter','latex')
            stackedplot(obj.T, obj.X(5:7,:)' )
            
            
            subplot(3,3,3);
            hold on
            grid on
            stairs(obj.T',[obj.U(1,1:end-1), obj.U(1,end-1)]','.-','Color',color)
            %xlabel('Step $k$','Interpreter','latex')
            subtitle('Control input ${u}$','Interpreter','latex');
            ylabel('${u}_1$','Interpreter','latex')
                        set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'fontsize',7);
            
            subplot(3,3,6);
            hold on
            grid on
            stairs(obj.T',[obj.U(2,1:end-1), obj.U(2,end-1)]','.-','Color',color)
            %xlabel('Step $k$','Interpreter','latex')
            set(gca,'fontsize',7);
            ylabel('${u}_2$','Interpreter','latex')
            
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            
            subplot(3,3,9);
            hold on
            grid on
            stairs(obj.T',[obj.U(3,1:end-1), obj.U(3,end-1)]','.-','Color',color)
            xlabel('Step $k$','Interpreter','latex')
            set(gca,'fontsize',7);
            ylabel('${u}_3$','Interpreter','latex')
            
        end
        function obj = plot_area(obj, X_low, X_up, U_low, U_up,color, fig_num,txt)
            
            figure(fig_num);
            subplot(6,2,2);
            hold on
            grid on
            fill_area(obj.T, X_up(5,:), X_low(5,:), color,txt);
            subplot(6,2,4);
            hold on
            grid on
            fill_area(obj.T, X_up(6,:), X_low(6,:), color,txt);
            
            subplot(6,2,6);
            hold on
            grid on
            fill_area(obj.T, X_up(7,:), X_low(7,:), color,txt);
            
            
            subplot(2,3,4);
            hold on
            grid on
            [x1, y1] = stairs(obj.T, [U_up(1,:), U_up(1,end)]);
            [x2, y2] = stairs(obj.T, [U_low(1,:), U_low(1,end)]);
            fill_area(x1', y1', y2', color,txt);
            
            
            subplot(2,3,5);
            hold on
            grid on
            [x1, y1] = stairs(obj.T, [U_up(2,:), U_up(2,end)]);
            [x2, y2] = stairs(obj.T, [U_low(2,:), U_low(2,end)]);
            fill_area(x1', y1', y2', color,txt);
            
            subplot(2,3,6);
            hold on
            grid on
            
            [x1, y1] = stairs(obj.T, [U_up(3,:), U_up(3,end)]);
            [x2, y2] = stairs(obj.T, [U_low(3,:), U_low(3,end)]);
            fill_area(x1', y1', y2', color,txt);
        end
        
        function obj = plot_tubes_phase(obj, delta ,color,txt,traj_noise,s)
            %M = [1,9;11,19];
            %s = subplot(2,10,M(fig_num,:));
            title(txt);
            set(gca,'fontsize',7);
            hold on;
            
            p = plot(obj.X(7,1:end)+delta(7,1:end),obj.X(6,1:end)+delta(6,1:end),'--','Color',color);
            plot(obj.X(7,1:end)+delta(7,1:end),obj.X(6,1:end)-delta(7+6,1:end),'--','Color',color);
            plot( obj.X(7,1:end)-delta(7+7,1:end),obj.X(6,1:end)+delta(6,1:end),'--','Color',color);
            plot( obj.X(7,1:end)-delta(7+7,1:end),obj.X(6,1:end)-delta(7+6,1:end),'--','Color',color);
            
            c = turbo(11);
            for i=obj.K:-1:1
                x = [obj.X(7,i)+delta(7,i), obj.X(7,i)+delta(7,i), obj.X(7,i)-delta(7+7,i), obj.X(7,i)-delta(7+7,i)];
                y = [obj.X(6,i)+delta(6,i), obj.X(6,i)-delta(7+6,i), obj.X(6,i)-delta(7+6,i) , obj.X(6,i)+delta(6,i)];
                s = fill(x,y, c(i,:),'DisplayName',txt);
                s.EdgeColor = color;
                alpha(s,.3);
            end
            rectangle('Position',[-obj.m.w_max -obj.m.w_max 2*obj.m.w_max 2*obj.m.w_max],'EdgeColor','k','LineStyle','-','Linewidth',1.5);
            axis square
            plot( obj.X(7,:),obj.X(6,:),'.-','Color',"#0072BD",'Linewidth',1.5);%blue
            plot( traj_noise.X(7,:),traj_noise.X(6,:),'.-','Color',"#D95319",'Linewidth',1.5);%red
            
        end
        
        function obj = plot_tubes(obj, delta ,color, fig_num,txt)
            
            figure(fig_num);
            subplot(3,3,[1,4]);
            subtitle('Quaternion ${q}$','Interpreter','latex')
            hold on
            grid on
            fill_area(obj.T,obj.X(1,:)+delta(1,:),obj.X(1,:)-delta(1+7,:), color,txt);
            fill_area(obj.T,obj.X(2,:)+delta(2,:),obj.X(2,:)-delta(2+7,:), color,txt);
            fill_area(obj.T,obj.X(3,:)+delta(3,:),obj.X(3,:)-delta(3+7,:), color,txt);
            fill_area(obj.T,obj.X(4,:)+delta(4,:),obj.X(4,:)-delta(4+7,:), color,txt);
            
            
            subplot(3,3,2);
            hold on
            grid on
            fill_area(obj.T, obj.X(5,:)+delta(5,:), obj.X(5,:)-delta(5+7,:), color,txt);
            subplot(3,3,5);
            hold on
            grid on
            fill_area(obj.T, obj.X(6,:)+delta(6,:), obj.X(6,:)-delta(6+7,:),color,txt);
            subplot(3,3,8);
            hold on
            grid on
            fill_area(obj.T, obj.X(7,:)+delta(7,:), obj.X(7,:)-delta(7+7,:), color,txt);
            
            subplot(3,3,3);
            hold on
            grid on
            [x1, y1] = stairs(obj.T, [obj.U(1,:)]+delta(14+1,:));
            [x2, y2] = stairs(obj.T, [obj.U(1,:)]-delta(14+4,:));
            x1 = x1(1:end-1);
            x2 = x2(1:end-1);
            y1 = y1(1:end-1);
            y2 = y2(1:end-1);
            fill_area(x1', y1', y2', color,txt);
            
            subplot(3,3,6);
            hold on
            grid on
            [x1, y1] = stairs(obj.T, [obj.U(2,:)]+delta(14+2,:));
            [x2, y2] = stairs(obj.T, [obj.U(2,:)]-delta(14+5,:));
            x1 = x1(1:end-1);
            x2 = x2(1:end-1);
            y1 = y1(1:end-1);
            y2 = y2(1:end-1);
            fill_area(x1', y1', y2', color,txt);
            
            subplot(3,3,9);
            hold on
            grid on
            [x1, y1] = stairs(obj.T, [obj.U(3,:)]+delta(14+3,:));
            [x2, y2] = stairs(obj.T, [obj.U(3,:)]-delta(14+6,:));
            x1 = x1(1:end-1);
            x2 = x2(1:end-1);
            y1 = y1(1:end-1);
            y2 = y2(1:end-1);
            fill_area(x1', y1', y2',color,txt);
        end
        
        
        function obj = plot_quat(obj,fig_num)
            
            figure(fig_num);
            hold on
            grid on
            for i = 1:obj.K
                eu_1(:,i) = rot_quat((obj.X(1:4,i)), [1;0;0]);
                eu_2(:,i) = rot_quat((obj.X(1:4,i)), [0;1;0]);
                eu_3(:,i) = rot_quat((obj.X(1:4,i)), [0;0;1]);
            end
            plot3(eu_1(1,:),eu_1(2,:),eu_1(3,:),'Linewidth',4,'color','r');
            plot3(eu_2(1,:),eu_2(2,:),eu_2(3,:),'Linewidth',4,'color','b');
            plot3(eu_3(1,:),eu_3(2,:),eu_3(3,:),'Linewidth',4,'color','g');
            
            xlabel('X-axis');
            ylabel('Y-axis');
            zlabel('Z-axis');
            [X,Y,Z] = sphere(100);
            surf(X,Y,Z,'FaceAlpha',0.5);
            axis equal;
        end
        
        
        
        function l = add_constraints(obj,model,fig_num)
            figure(fig_num);
            subplot(3,3,2);
            hold on;
            yline(model.w_max,'k-', 'Linewidth',1.5,'DisplayName','Constraints');
            yline(-model.w_max,'k-', 'Linewidth',1.5);
            %l = legend;
            h=get(gca,'Children');
            legendstr={'Constraints','Non-robust.','LTV','Uncertain sys','Nonlinear reach. set','Nominal'};
            l = legend(h([2 3 5 4 6 7]),legendstr{[1 2 3 4 5 6]},'Interpreter','latex');
            l.Position = [0.1450 0.0651 0.1629 0.2776];
            l.Box = 'off'; % Remove the legend box
            l.FontSize = 6; % Adjust the font size as needed
            l.NumColumns = 2; % Set the number of legend columns to 2
            xylim = 1.5;
            xylim_input = 1.05;
            
            ylim([-xylim*model.w_max, xylim*model.w_max])
            set(gca,'FontName','Times')
            
            subplot(3,3,5);
            hold on;
            yline(model.w_max,'k-', 'Linewidth',1.5);
            yline(-model.w_max,'k-', 'Linewidth',1.5);
            ylim([-xylim*model.w_max, xylim*model.w_max])
            set(gca,'FontName','Times')
            
            subplot(3,3,8);
            hold on;
            yline(model.w_max,'k-', 'Linewidth',1.5);
            yline(-model.w_max,'k-', 'Linewidth',1.5);
            ylim([-xylim*model.w_max, xylim*model.w_max])
            set(gca,'FontName','Times')
            
            
            subplot(3,3,3);
            hold on
            yline(model.T_max,'k-', 'Linewidth',1.5);
            yline(-model.T_max,'k-', 'Linewidth',1.5);
            ylim([-xylim_input*model.T_max, xylim_input*model.T_max])
            set(gca,'FontName','Times')
            
            
            subplot(3,3,6);
            hold on
            yline(model.T_max,'k-', 'Linewidth',1.5);
            yline(-model.T_max,'k-', 'Linewidth',1.5);
            ylim([-xylim_input*model.T_max, xylim_input*model.T_max])
            set(gca,'FontName','Times')
            
            
            subplot(3,3,9);
            hold on
            yline(model.T_max,'k-', 'Linewidth',1.5);
            yline(-model.T_max,'k-', 'Linewidth',1.5);
            ylim([-xylim_input*model.T_max, xylim_input*model.T_max])
            set(gca,'FontName','Times')
            
        end
        
    end
end

