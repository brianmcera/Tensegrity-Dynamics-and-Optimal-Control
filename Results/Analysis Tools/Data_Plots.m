%% Plotting and Data Analysis #########################################
%######################################################################
%t_time: last timestep to simulate
%t_start: first timestep to simulate

%% COM trajectory
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last')

t_start = 1;

drawBasePolygons = 1;

baseFloor = min(X_record.p(3:3:end,1))+1e-2;

% des_direction = simulationParameters_record.costArgs.target;
dt = simulationParameters_record.timestep;

COM_x = mean(X_record.p(1:3:end,:));
COM_y = mean(X_record.p(2:3:end,:));

figure(1)
plot(COM_x(t_start:t_time),COM_y(t_start:t_time),'b-','linewidth',4)
hold on
plot(COM_x(t_start),COM_y(t_start),'rp','linewidth',4)
speed = norm([COM_x(t_time),COM_y(t_time)]-...
    [COM_x(t_start),COM_y(t_start)])/(dt*(t_time-t_start));
title(['COM of Robot over Time, Avg. Velocity = ',num2str(speed),' m/s'])
axis('equal')
xlabel('X [m]')
ylabel('Y [m]')

if(drawBasePolygons)
    for i = 1:20:t_time
        nodesXYZ = [X_record.p(1:3:end,i),...,
            X_record.p(2:3:end,i),...
            X_record.p(3:3:end,i)];
        nodesXYZ = nodesXYZ(X_record.p(3:3:end,i)<=baseFloor,1:2);
        patch('XData',nodesXYZ(:,1),'YData',nodesXYZ(:,2),...
            'EdgeColor','b','FaceColor','none','LineWidth',0.5,'LineStyle','--');
        pause(1e-3)
    end
end

%custom code
% scatter([0 0 3 3],[0 -3 -3 0],'gp','linewidth',8)
% scatter([0],[0],'rp','linewidth',15)
% rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'FaceColor',[0.5 1 0.5 0.25])
% rectangle('Position',[-0.5 -3.5 1 1],'Curvature',[1 1],'FaceColor',[0.5 1 0.5 0.25])
% rectangle('Position',[2.5 -3.5 1 1],'Curvature',[1 1],'FaceColor',[0.5 1 0.5 0.25])
% rectangle('Position',[2.5 -0.5 1 1],'Curvature',[1 1],'FaceColor',[0.5 1 0.5 0.25])

set(gca,'fontsize',20,'fontname','times new roman')
grid on
    
%% 3D COM trajectory
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last')

t_start = 1;

% des_direction = simulationParameters_record.costArgs.target;
dt = simulationParameters_record.timestep;

COM_x = mean(X_record.p(1:3:end,:));
COM_y = mean(X_record.p(2:3:end,:));
COM_z = mean(X_record.p(3:3:end,:));

figure(1)
plot3(COM_x(t_start:t_time),COM_y(t_start:t_time),COM_z(t_start:t_time),...
    'linewidth',2)
hold on
speed = norm([COM_x(t_time),COM_y(t_time)]-...
    [COM_x(t_start),COM_y(t_start)])/(dt*(t_time-t_start));
title(['COM of Robot over Time, Avg. Velocity = ',...
    num2str(speed)])
% disp(['Distance along Desired Direction:',num2str(([COM_x(t_time),...
%     COM_y(t_time)]-[COM_x(t_start),COM_y(t_start)])*des_direction(1:2)')])
xlabel('x')
ylabel('y')
zlabel('z')
axis('equal')

%% Dynamic movement WITHOUT ghost trail image
h = figure(1);
clf
title('Press Any Key to Start')
pause()

isolateNodes = false; %only plot nodes
plot_openLoopTraj = false; %MPC OpenLoop Trajectories for each node
highlightSupportNodes = false;
plotDirection = true;
plotCOM = true;
recordVideo = false;

if(recordVideo)
    v = VideoWriter(['Rolling_CentralPayload_',datestr(now,'yy-mm-dd_HH_MM_SS'),...
        '.avi']);
    v.FrameRate = 33;  % Default 30
    v.Quality = 75;    % Default 75
    open(v);
end


rotateCamera = false; %continuous rotation camera
el = 10; %elevation camera angle
az = 0; %azimuthal camera angle

t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');

omega = simulationParameters_record.omega;
constraints = [];

for i = t_start:10:t_time
    %pause()
    if(~isolateNodes)
        Xbar.p = X_record.p(:,i);
        Xbar.pDOT = X_record.pDOT(:,i);
        structurePlot(Xbar,omega,constraints,[az,el],1,1,0)
        camlight

        if(rotateCamera)
            az = az+2;
        end
        title(['time = ',num2str(i*dt)])

    else
        clf %for the node/direction plot
    end
    
    %plot floor
    CoM = [mean(X_record.p(1:3:end,i));
        mean(X_record.p(2:3:end,i));
        mean(X_record.p(3:3:end,i))];
    fill3([CoM(1)+1;CoM(1)+1;CoM(1)-1;CoM(1)-1],...
        [CoM(2)+1;CoM(2)-1;CoM(2)-1;CoM(2)+1],...
        ones(4,1)*baseFloor,[0.3 0.3 0.3],'facealpha',0.5)
    
    if(plot_openLoopTraj)        
        %plot open loop MPC trajectories for each node
        hold on
        p_OL = X_openLoop_record{i}.p(:,1:end-1);
        for k = 1:size(p_OL,1)/3
            plot3(p_OL(k*3-2,:),p_OL(k*3-1,:),p_OL(k*3,:),'y-','linewidth',3);
            plot3(p_OL(k*3-2,end),p_OL(k*3-1,end),p_OL(k*3,end),'rp','linewidth',1);
        end
        hold on
        
        lead_direction = [controllerOutputArgs_record{i}.desDirection ];
        %plot desired roll direction for each node
        Centroid_pos = [mean(p_OL(1:3:end,1));...
            mean(p_OL(2:3:end,1));...
            mean(p_OL(3:3:end,1))];
        cross_axis = cross(lead_direction,[0 0 1]);
        cross_axis = cross_axis/norm(cross_axis);
        for node = 1:length(p_OL(:,1))/3
            %
            %cross product velocity reward
            r_vec = p_OL(node*3-2:node*3,1)-Centroid_pos;
            rolling_dir = cross(r_vec,cross_axis);
            rolling_dir = rolling_dir*0.5; %scaling for graphics
            %rolling_dir = [0 0 1]*norm(rolling_dir)+rolling_dir;
            %rolling_dir = rolling_dir/norm(rolling_dir);
            quiver3(p_OL(node*3-2,1),p_OL(node*3-1,1),p_OL(node*3,1),...
                rolling_dir(1),rolling_dir(2),rolling_dir(3),'b','linewidth',1);
        end
        
        if(isolateNodes)
            view(0,0)
        end
    end
    baseFloor = min(X_record.p(3:3:end,1))+5e-3;
    if(highlightSupportNodes)
        for node = 1:length(X_record.p(:,1))/3
            if(X_record.p(node*3,i)<=baseFloor)
                plot3(X_record.p(node*3-2,i),X_record.p(node*3-1,i),...
                    X_record.p(node*3,i),'bsq','linewidth',3);
            end
        end
    end
    
    %highlight one node
    hNode = 7;
    plot3(X_record.p(hNode*3-2,i),X_record.p(hNode*3-1,i),...
        X_record.p(hNode*3,i),'rp','linewidth',3);
    
    %plot desired direction
    if(plotDirection)
        Centroid_pos = [mean(X_record.p(1:3:end,i));...
            mean(X_record.p(2:3:end,i));...
            mean(X_record.p(3:3:end,i))];
        lead_direction = [controllerOutputArgs_record{i}.desDirection];
        quiver3(Centroid_pos(1),Centroid_pos(2),Centroid_pos(3),...
            lead_direction(1)/5,lead_direction(2)/5,lead_direction(3)/5,'b','linewidth',1);
    end
    
    %plot desired direction
    if(plotCOM)
        Centroid_pos = [mean(X_record.p(1:3:end,i));...
            mean(X_record.p(2:3:end,i));...
            mean(X_record.p(3:3:end,i))];
        plot3([Centroid_pos(1) Centroid_pos(1)],...
            [Centroid_pos(2) Centroid_pos(2)],...
            [Centroid_pos(3) baseFloor],'r--p','linewidth',5);
    end
    
    zlim([0 1])
    pause(dt/10)
    
    if(recordVideo)
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end

if(recordVideo)
    close(v)
end
   
%% Dynamic movement, ODE Matrix Input 
h = figure(1);
clf
title('Press Any Key to Start')
pause()

isolateNodes = false; %only plot nodes
plot_openLoopTraj = false; %MPC OpenLoop Trajectories for each node
highlightSupportNodes = true;
plotDirection = true;
plotCOM = true;
recordVideo = false;

rotateCamera = false; %continuous rotation camera
el = 15; %elevation camera angle
az = 0; %azimuthal camera angle
constraints = [];
omega = I_six_bar_model()
xTrans = XOUT';
for idx = 1:size(xTrans,2)
    %pause()
    cla %for the node/direction plot

%     CoM = [mean(xTrans(1:3:36,idx));
%         mean(xTrans(2:3:36,idx));
%         mean(xTrans(3:3:36,idx))];
    
    Xbar.p = xTrans(1:36,idx);
    structurePlot(Xbar,omega,[],[az,el],1,1,1)

    pause(1e-3)
    
end

%% Dynamics Movement, State Estimates
h = figure(1);
clf
title('Press Any Key to Start')
pause()

isolateNodes = false; %only plot,  nodes
plot_openLoopTraj = false; %MPC OpenLoop Trajectories for each node
highlightSupportNodes = true;
plotDirection = true;
plotCOM = true;
plotUnc = true;


rotateCamera = false; %continuous rotation camera
el = 15; %elevation camera angle
az = 0; %azimuthal camera angle

t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((Xhat_record.p(1,:)),1,'last');

omega = simulationParameters_record.omega;
constraints = [];

for i = t_start:1:t_time
    %pause()
    
    subplot(2,2,2)
    cla
    bar([sqrt(Pm_record(1:36,i));abs([Xhat_record.p(:,i)-X_record.p(:,i)])])
    title('Position State Estimates')
    %hold on
%     bar(abs([Xhat_record.p(:,i)-X_record.p(:,i)]))
    subplot(2,2,4)
    cla
    bar([sqrt(Pm_record(73:96,i));abs([Xhat_record.RL(:,i)-X_record.RL(:,i)])])
    title('Cable State Estimates')
    if(~isolateNodes)
        Xbar.p = Xhat_record.p(:,i);
        Xbar.pDOT = Xhat_record.pDOT(:,i);
        subplot(2,2,[1 3])
        %cla
        structurePlot(Xbar,omega,constraints,[az,el],1,1,1)
        camlight

        if(rotateCamera)
            az = az+2;
        end
        title(['time = ',num2str(i*dt)])

    else
        clf %for the node/direction plot
    end
    
    %plot floor
    CoM = [mean(Xhat_record.p(1:3:end,i));
        mean(Xhat_record.p(2:3:end,i));
        mean(Xhat_record.p(3:3:end,i))];
    fill3([CoM(1)+1;CoM(1)+1;CoM(1)-1;CoM(1)-1],...
        [CoM(2)+1;CoM(2)-1;CoM(2)-1;CoM(2)+1],...
        ones(4,1)*baseFloor,[0.3 0.3 0.3],'facealpha',0.5)
    
    if(plot_openLoopTraj)        
        %plot open loop MPC trajectories for each node
        hold on
        p_OL = X_openLoop_record{i+1}.p;
        for k = 1:12
            plot3(p_OL(k*3-2,:),p_OL(k*3-1,:),p_OL(k*3,:),'y-','linewidth',3);
            plot3(p_OL(k*3-2,end),p_OL(k*3-1,end),p_OL(k*3,end),'rp','linewidth',1);
        end
        hold on
        
        lead_direction = [controllerOutputArgs_record{i}.desDirection ];
        %plot desired roll direction for each node
        Centroid_pos = [mean(p_OL(1:3:end,1));...
            mean(p_OL(2:3:end,1));...
            mean(p_OL(3:3:end,1))];
        cross_axis = cross(lead_direction,[0 0 1]);
        cross_axis = cross_axis/norm(cross_axis);
        for node = 1:length(p_OL(:,1))/3
            %
            %cross product velocity reward
            r_vec = p_OL(node*3-2:node*3,1)-Centroid_pos;
            rolling_dir = cross(r_vec,cross_axis);
            rolling_dir = rolling_dir*0.5; %scaling for graphics
            %rolling_dir = [0 0 1]*norm(rolling_dir)+rolling_dir;
            %rolling_dir = rolling_dir/norm(rolling_dir);
            quiver3(p_OL(node*3-2,1),p_OL(node*3-1,1),p_OL(node*3,1),...
                rolling_dir(1),rolling_dir(2),rolling_dir(3),'b','linewidth',1);
        end
        
        if(isolateNodes)
            view(0,0)
        end
    end
    baseFloor = -0.4;%min(min(Xhat_record.p(3:3:end,1:t_time)))+5e-3;
    if(highlightSupportNodes)
        for node = 1:length(Xhat_record.p(:,1))/3
            if(Xhat_record.p(node*3,i)<=baseFloor)
                plot3(Xhat_record.p(node*3-2,i),Xhat_record.p(node*3-1,i),...
                    Xhat_record.p(node*3,i),'bsq','linewidth',3);
            end
        end
    end
    
    %visualize nodal position uncertainty
    if(plotUnc)
        R = abs(simulationParameters_record.omega.R);
        for j = 1:size(R,1)
            idx = find(R(j,:));
            for k = 1:numel(idx)
            [X,Y,Z] = ellipsoid(Xhat_record.p(idx(k)*3-2,i),...
                Xhat_record.p(idx(k)*3-1,i),Xhat_record.p(idx(k)*3,i),...
                sqrt(Pm_record(idx(k)*3-2,i)),sqrt(Pm_record(idx(k)*3-1,i)),...
                sqrt(Pm_record(idx(k)*3,i)));
            surf(X,Y,Z,'FaceAlpha',0.2);
            end
        end
    end
    
    %highlight one node
    hNode = 7;
    plot3(Xhat_record.p(hNode*3-2,i),Xhat_record.p(hNode*3-1,i),...
        Xhat_record.p(hNode*3,i),'rp','linewidth',3);
    
    %plot desired direction
    if(plotDirection)
        Centroid_pos = [mean(Xhat_record.p(1:3:end,i));...
            mean(Xhat_record.p(2:3:end,i));...
            mean(Xhat_record.p(3:3:end,i))];
        lead_direction = [controllerOutputArgs_record{i}.desDirection];
        quiver3(Centroid_pos(1),Centroid_pos(2),Centroid_pos(3),...
            lead_direction(1)/5,lead_direction(2)/5,lead_direction(3)/5,'b','linewidth',5);
    end
    
    %plot desired direction
    if(plotCOM)
        Centroid_pos = [mean(Xhat_record.p(1:3:end,i));...
            mean(Xhat_record.p(2:3:end,i));...
            mean(Xhat_record.p(3:3:end,i))];
        plot3([Centroid_pos(1) Centroid_pos(1)],...
            [Centroid_pos(2) Centroid_pos(2)],...
            [Centroid_pos(3) baseFloor],'r--p','linewidth',5);
    end
    
    zlim([baseFloor-0.01 .4])
    pause(dt/10)
end 
 
%% Dynamics Movement, True State AND State Estimates
h = figure(1);
clf
title('Press Any Key to Start')
pause()

isolateNodes = false; %only plot nodes
plot_openLoopTraj = false; %MPC OpenLoop Trajectories for each node
highlightSupportNodes = true;
plotDirection = true;
plotCOM = true;
plotUnc = true;


rotateCamera = false; %continuous rotation camera
el = 0; %elevation camera angle
az = 0; %azimuthal camera angle

t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((Xhat_record.p(1,:)),1,'last');

omega = simulationParameters_record.omega;
constraints = [];

for i = t_start:5:t_time
    %pause()
    if(~isolateNodes)
        Xhat.p = Xhat_record.p(:,i);
        Xhat.pDOT = Xhat_record.pDOT(:,i);
        Xbar.p = X_record.p(:,i);
        Xbar.pDOT = X_record.pDOT(:,i);
        cla
        structurePlot(Xhat,omega,constraints,[az,el],0,1,1)
        structurePlot(Xbar,omega,constraints,[az,el],0,0,0)

        drawnow
        camlight

        if(rotateCamera)
            az = az+2;
        end
        title(['time = ',num2str(i*dt)])

    else
        clf %for the node/direction plot
    end
    
    %plot floor
    CoM = [mean(X_record.p(1:3:end,i));
        mean(X_record.p(2:3:end,i));
        mean(X_record.p(3:3:end,i))];
    fill3([CoM(1)+1;CoM(1)+1;CoM(1)-1;CoM(1)-1],...
        [CoM(2)+1;CoM(2)-1;CoM(2)-1;CoM(2)+1],...
        ones(4,1)*baseFloor,[0.3 0.3 0.3],'facealpha',0.5)
    
    if(plot_openLoopTraj)        
        %plot open loop MPC trajectories for each node
        hold on
        p_OL = X_openLoop_record{i+1}.p;
        for k = 1:12
            plot3(p_OL(k*3-2,:),p_OL(k*3-1,:),p_OL(k*3,:),'y-','linewidth',3);
            plot3(p_OL(k*3-2,end),p_OL(k*3-1,end),p_OL(k*3,end),'rp','linewidth',1);
        end
        hold on
        
        lead_direction = [controllerOutputArgs_record{i}.desDirection ];
        %plot desired roll direction for each node
        Centroid_pos = [mean(p_OL(1:3:end,1));...
            mean(p_OL(2:3:end,1));...
            mean(p_OL(3:3:end,1))];
        cross_axis = cross(lead_direction,[0 0 1]);
        cross_axis = cross_axis/norm(cross_axis);
        for node = 1:length(p_OL(:,1))/3
            %
            %cross product velocity reward
            r_vec = p_OL(node*3-2:node*3,1)-Centroid_pos;
            rolling_dir = cross(r_vec,cross_axis);
            rolling_dir = rolling_dir*0.5; %scaling for graphics
            %rolling_dir = [0 0 1]*norm(rolling_dir)+rolling_dir;
            %rolling_dir = rolling_dir/norm(rolling_dir);
            quiver3(p_OL(node*3-2,1),p_OL(node*3-1,1),p_OL(node*3,1),...
                rolling_dir(1),rolling_dir(2),rolling_dir(3),'b','linewidth',1);
        end
        
        if(isolateNodes)
            view(0,0)
        end
    end
    baseFloor = min(X_record.p(3:3:end,1))+5e-3;
    if(highlightSupportNodes)
        for node = 1:length(Xhat_record.p(:,1))/3
            if(Xhat_record.p(node*3,i)<=baseFloor)
                plot3(Xhat_record.p(node*3-2,i),Xhat_record.p(node*3-1,i),...
                    Xhat_record.p(node*3,i),'bsq','linewidth',3);
            end
        end
    end
    
    %visualize nodal position uncertainty
    if(plotUnc)
        R = abs(simulationParameters_record.omega.R);
        for j = 1:size(R,1)
            idx = find(R(j,:));
            for k = 1:numel(idx)
            [X,Y,Z] = ellipsoid(Xhat_record.p(idx(k)*3-2,i),...
                Xhat_record.p(idx(k)*3-1,i),Xhat_record.p(idx(k)*3,i),...
                3*sqrt(Pm_record(idx(k)*3-2,i)),3*sqrt(Pm_record(idx(k)*3-1,i)),...
                3*sqrt(Pm_record(idx(k)*3,i)));
            surf(X,Y,Z,'FaceAlpha',0.2);
            end
        end
    end
    
    %highlight one node
    hNode = 7;
    plot3(Xhat_record.p(hNode*3-2,i),Xhat_record.p(hNode*3-1,i),...
        Xhat_record.p(hNode*3,i),'rp','linewidth',3);
    
    %plot desired direction
    if(plotDirection)
        Centroid_pos = [mean(Xhat_record.p(1:3:end,i));...
            mean(Xhat_record.p(2:3:end,i));...
            mean(Xhat_record.p(3:3:end,i))];
        lead_direction = [controllerOutputArgs_record{i}.desDirection];
        quiver3(Centroid_pos(1),Centroid_pos(2),Centroid_pos(3),...
            lead_direction(1)/5,lead_direction(2)/5,lead_direction(3)/5,'b','linewidth',5);
    end
    
    %plot desired direction
    if(plotCOM)
        Centroid_pos = [mean(Xhat_record.p(1:3:end,i));...
            mean(Xhat_record.p(2:3:end,i));...
            mean(Xhat_record.p(3:3:end,i))];
        plot3([Centroid_pos(1) Centroid_pos(1)],...
            [Centroid_pos(2) Centroid_pos(2)],...
            [Centroid_pos(3) baseFloor],'r--p','linewidth',5);
    end
    
    %zlim([baseFloor-0.01 .4])
    pause(dt/10)
    pDiff = norm(Xhat_record.p(:,i)-X_record.p(:,i))
    RLDiff = norm(Xhat_record.RL(:,i)-X_record.RL(:,i),1)/24
end 

%% Dynamic movement WITH ghost trail image
h = figure(1);
clf
title('Press Any Key to Start')
pause()

isolateNodes = false; %only plot nodes
plot_openLoopTraj = false; %MPC OpenLoop Trajectories for each node
highlightSupportNodes = true;
plotDirection = false;

rotateCamera = false; %continuous rotation camera
el = 20; %elevation camera angle
az = 45; %azimuthal camera angle

t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');

omega = simulationParameters_record.omega;
constraints = [];

firstIter = 1;
batches = 3;
for i = t_start+1:10:t_time

    if(~isolateNodes)
        Xbar.p = X_record.p(:,i);
        Xbar.pDOT = X_record.pDOT(:,i);
        structurePlot(Xbar,omega,constraints,[az,el],0,1)
        %         h2 = get(gca,'Children');
        h2 = findobj(gca,'Type','Line');
        if(firstIter == 1 && numel(h2)~=0)
            numLines = numel(h2);
            firstIter = 0;
        elseif(numel(h2)~=0)
            for line = 1:numel(h2)
                batchNum = ceil(line/numLines)
                if(batchNum~=1)
                    h2(line).Color(1:3) =[0 0 1];
                    h2(line).Color(4) =(1/(exp(batchNum*0.7)));
                end
            end
        end
        if(numel(h2)>(batches-1)*numLines)
            delete(h2(numel(h2)-numLines+1:end))
            %TODO: CURRENTLY EXPERIENCING BUG WHERE H2 ALSO KEEPS TRACK OF
            %BLUE SQUARE MARKERS. THUS, THE SIZE OF H2 IS NOT CONSTANT OVER
            %TIME AND IS DEPENDENT ON THE CONTACT CONDITIONS
        end
        delete(findobj(gca,'Type','Scatter'))
                
        if(rotateCamera)
            az = az+2;
        end
        title(['time = ',num2str(i*dt)])

    else
        clf %for the node/direction plot
    end
    
    if(plot_openLoopTraj)        
        %plot open loop MPC trajectories for each node
        hold on
        p_OL = X_openLoop_record{i+1}.p;
        for k = 1:12
            plot3(p_OL(k*3-2,:),p_OL(k*3-1,:),p_OL(k*3,:),'y-','linewidth',3);
            plot3(p_OL(k*3-2,end),p_OL(k*3-1,end),p_OL(k*3,end),'rp','linewidth',1);
        end
        hold on
        
        lead_direction = [controllerOutputArgs_record{i}.desDirection ];
        %plot desired roll direction for each node
        Centroid_pos = [mean(p_OL(1:3:end,1));...
            mean(p_OL(2:3:end,1));...
            mean(p_OL(3:3:end,1))];
        cross_axis = cross(lead_direction,[0 0 1]);
        cross_axis = cross_axis/norm(cross_axis);
        for node = 1:length(p_OL(:,1))/3
            %
            %cross product velocity reward
            r_vec = p_OL(node*3-2:node*3,1)-Centroid_pos;
            rolling_dir = cross(r_vec,cross_axis);
            rolling_dir = rolling_dir*0.5; %scaling for graphics
            %rolling_dir = [0 0 1]*norm(rolling_dir)+rolling_dir;
            %rolling_dir = rolling_dir/norm(rolling_dir);
            quiver3(p_OL(node*3-2,1),p_OL(node*3-1,1),p_OL(node*3,1),...
                rolling_dir(1),rolling_dir(2),rolling_dir(3),'b','linewidth',1);
        end
        
        if(isolateNodes)
            view(0,0)
        end
    end
    baseFloor = min(X_record.p(3:3:end,1))+5e-3;
    if(highlightSupportNodes)
        for node = 1:length(X_openLoop_record{1}.p)/3
            if(X_record.p(node*3,i)<=baseFloor)
                scatter3(X_record.p(node*3-2,i),X_record.p(node*3-1,i),...
                    X_record.p(node*3,i),'bsq','linewidth',3);
            end
        end
    end
    
    pause(dt/10)
end 
 
%% Cable Lengths
%only look at time steps up to when simulation was aborted
figure(2)
dt = simulationParameters_record.timestep;
%time up to   last simulated time
t_time = find((X_record.p(1,:)),1,'last')
t_gap = 1;
%color_list = hsv(32);
ax1 = subplot(2,1,1);

tFin = t_time;
t_start = tFin;
%pause()

for t_time = t_start:t_gap:tFin
    clf
    for i = 1:size(X_record.RL,1)
        subplot(2,1,1);
        plot((0:t_time-1)*dt,X_record.RL(i,1:t_time)*100,'-','linewidth',2)%,'col',color_list(i,:))
        hold on
    end
    pause(1e-4)
    t_time
end
colormap(hot)
%title('Absolute Cable Restlengths','fontsize',20)
grid on
xlabel('Time [s]')
ylabel('Cable Rest Length [cm]')
set(gca,'fontsize',30,'fontname','times new roman')


% ax2 = subplot(2,1,2); title('Restlength Changes from Pretension [m]');
% for i = 1:size(U_record.RL,1)
%     plot((0:t_time-1)*dt,(U_record.RL(i,1:t_time)-U_record.RL(i,1)))%,'col',color_list(i,:))
%     hold on
% end
% linkaxes([ax1,ax2],'x')
% title('Restlength Changes from Initial Condition','fontsize',20)
% ylabel('Relative Restlengths [m]','fontsize',20)
% xlabel('Time [s]','fontsize',20)
% grid on
%  

%% Cable Velocity Inputs
%only look at time steps up to when simulation was aborted
figure(2)
dt = simulationParameters_record.timestep;
%time up to   last simulated time
t_time = find((X_record.p(1,:)),1,'last')
t_gap = 1;
%color_list = hsv(32);
ax1 = subplot(2,1,1);

tFin = t_time;
t_start = tFin;
%pause()

for t_time = t_start:t_gap:tFin
    clf
    for i = 1:size(X_record.RL,1)
        subplot(2,1,1);
        plot((0:t_time-1)*dt,U_record.RLdot(i,1:t_time)*100,'-','linewidth',2)%,'col',color_list(i,:))
        hold on
    end
    pause(1e-4)
    t_time
end
colormap(hot)
%title('Absolute Cable Restlengths','fontsize',20)
grid on
xlabel('Time [s]')
ylabel('Cable Velocity [cm/s]')
set(gca,'fontsize',30,'fontname','times new roman')

%% Rod Lengths
%only look at time steps up to when simulation was aborted
figure(3)
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');
%color_list = hsv(32);
for i = 1:size(X_record.L,1)
    plot((0:t_time-1),X_record.L(i,1:t_time)*100)%,'col',color_list(i,:))
    hold on
end
title('Rod Lengths','fontsize',20)
ylabel('Relative Restlengths [cm]','fontsize',20)
xlabel('Time [s]','fontsize',20)

%% Cable Restlengths in bar plot
%only look at time steps up to when simulation was aborted
figure(4)
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');

clf
color_list = prism(32);
set(gca,'Color','k')
upper_limit = max(max(X_record.RL(:,1:t_time)));
lower_limit = min(min(X_record.RL(:,1:t_time)));
ax1 = subplot(2,1,1);
bar(X_record.RL(:,1),'b')
for j = 1:5:t_time
%     pause()
    bar(X_record.RL(:,j),'b')
    ylim([lower_limit upper_limit])
    title(['t = ',num2str(dt*j)])
    pause(1e-2)
end

%% Rod Lengths in bar plot
%only look at time steps up to when simulation was aborted
figure(5)
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');

clf
color_list = prism(32);
set(gca,'Color','k')
upper_limit = max(max(X_record.L(:,1:t_time)));
lower_limit = 0;
ax1 = subplot(2,1,1);
title('Rod Lengths over Time','fontsize',20)
for j = 1:1:t_time
    bar(X_record.L(:,j))
    ylim([lower_limit upper_limit])
    title(['t = ',num2str(dt*j)])
    pause(1e-2)
end
title('Rod Lengths over Time','fontsize',20)
ylabel('Rod Lengths [cm]','fontsize',20)
xlabel('Time [s]','fontsize',20)

%% Cable Tensions
figure(3)
subplot(2,1,1)
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');

C = simulationParameters_record.omega.C;
tensions = zeros(size(C,1),t_time);
stiffness = simulationParameters_record.omega.cables.stiffness;
for j = 1:t_time
    gamma = X_record.p(:,j);
    sep_dist = (kron(C,[1 0 0])*gamma).^2+...
        (kron(C,[0 1 0])*gamma).^2+...
        (kron(C,[0 0 1])*gamma).^2;
    sep_dist = sqrt(sep_dist);
    tensions(:,j) = stiffness.*(sep_dist-X_record.RL(:,j));
end
tensions(tensions<0)=0;
for cable = 1:size(C,1)%[1,4,21,9,12,18,17,16,10,22,8,2]%1:size(C,1)
    plot((1:t_time)*dt,tensions(cable,1:t_time))
    hold on
end
title('Cable Tensions Over Time')
xlabel('time [s]')
ylabel('tension [N]')
grid on

subplot(2,1,2)
plot((1:t_time)*dt,sum(tensions(:,1:t_time)))
grid on

%% Cable Tensions II
fig = figure(3)
left_color = [0.5 0.7 0.5];
right_color = [1 0 0];

colormap(pink)

ax1 = subplot(2,1,1)
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');

yyaxis left
C = simulationParameters_record.omega.C;
tensions = zeros(size(C,1),t_time);
stiffness = simulationParameters_record.omega.cables.stiffness;
for j = 1:t_time
    gamma = X_record.p(:,j);
    sep_dist = (kron(C,[1 0 0])*gamma).^2+...
        (kron(C,[0 1 0])*gamma).^2+...
        (kron(C,[0 0 1])*gamma).^2;
    sep_dist = sqrt(sep_dist);
    tensions(:,j) = stiffness.*(sep_dist-U_record.RL(:,j));
end
tensions(tensions<0)=0;
for cable = 1:size(C,1)
    plot((1:t_time)*dt,tensions(cable,1:t_time),'-')
    hold on
end
pairs = simulationParameters_record.omega.cables.paired;
plot((1:t_time)*dt,repmat(mean(mean(tensions(pairs,1:t_time))),1,t_time),'g-','linewidth',4)
xlabel('Time [s]')
ylabel('Cable Tension [N]')

colormap(winter)




yyaxis right
%plot((1:t_time)*dt,sum(tensions(:,1:t_time)),'k','linewidth',3)
h = area((1:t_time)*dt,sum(tensions(:,1:t_time)),(sum(tensions(:,1))))
h.FaceColor = [1,0,0];
h.FaceAlpha = 0.4;
h.EdgeColor = [1,0,0];
h.LineWidth = 2;
ylabel('Total Tension [N]')
%set(gca,'fontsize',36,'fontname','times new roman')
grid on

subplot(2,1,2)
pairs = simulationParameters_record.omega.cables.paired;
pair_tensions = zeros(size(pairs,1),size(tensions,2));
for i = 1:size(pairs,1)
    i
    pairs(i,:)
    pair_tensions(i,:) = abs(tensions(pairs(i,1),:)-tensions(pairs(i,2),:));
end
for pair = 1:size(pair_tensions,1)
    plot((1:t_time)*dt,pair_tensions(pair,1:t_time),'-')
    hold on
end
plot((1:t_time)*dt,repmat(mean(mean(pair_tensions(pair,1:t_time))),1,t_time),'g-','linewidth',4)

%% Power Draw
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last');

C = simulationParameters_record.omega.C;
tensions = zeros(size(C,1),t_time);
stiffness = simulationParameters_record.omega.cables.stiffness;
for j = 1:t_time
    gamma = X_record.p(:,j);
    sep_dist = (kron(C,[1 0 0])*gamma).^2+...
        (kron(C,[0 1 0])*gamma).^2+...
        (kron(C,[0 0 1])*gamma).^2;
    sep_dist = sqrt(sep_dist);
    tensions(:,j) = stiffness.*(sep_dist-U_record.RL(:,j));
end
tensions(tensions<0)=0;

cableDiff = diff(U_record.RL(:,1:t_time),1,2); %cable actuation

%handle paired-cables 
pairs = simulationParameters_record.omega.cables.paired;
for i = 1:size(pairs,1)
    tensions(pairs(i,1),:) = tensions(pairs(i,1),:)-tensions(pairs(i,2),:);
    %only count energy once for each pair
    tensions(pairs(i,2),:) = zeros(size(tensions(pairs(i,2),:)));
    cableDiff(pairs(i,2),:) = zeros(size(cableDiff(pairs(i,2),:)));
end
    
    
work = tensions(:,1:end-1).*-cableDiff; %force*dist (only when cables are retracted)
work(work<0) = 0;

workTime = sum(work);

total_work = sum(sum(work));

COM_x = mean(X_record.p(1:3:end,:));
COM_y = mean(X_record.p(2:3:end,:));

dist = norm([COM_x(t_time),COM_y(t_time)]-...
    [COM_x(t_start),COM_y(t_start)]);

motor_eff = 0.70;
COT = total_work/(.250*6)/dist/motor_eff

%% Cable Open Loop Trajectories
%only look at time steps up to when simulation was aborted
figure(4)
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last')

cable_nums = 25;
gap = 10;

for idx = 1:length(cable_nums)
    cable_num = cable_nums(idx);
    %plot actual trajectory
    plot((0:t_time-1)*dt,U_record.RL(cable_num,1:t_time),'linewidth',3)
    hold on
    
    %plot Open Loop trajectories from MPC
    title('Absolute Cable Restlengths','fontsize',20)
    
    N = simulationParameters_record.controllerHorizon;
    for t=1:gap:t_time
        plot((t-2:t-2+N)*dt,(U_openLoop_record{t}.RL(cable_num,1:N+1)),'k:')
        plot((t-2)*dt,(U_openLoop_record{t}.RL(cable_num,1)),'rp')
        hold on
    end
end

title('Restlength Changes from Initial Condition','fontsize',20)
ylabel('Relative Restlengths [m]','fontsize',20)
xlabel('Time [s]','fontsize',20)

%% Cost Results
figure(4)
t_start = 1;
dt = simulationParameters_record.timestep;
%time up to last simulated time
t_time = find((X_record.p(1,:)),1,'last')-1

costTrimmed = cost_record(1:t_time);
costSTD = std(costTrimmed);
numAnomalies = sum(abs(costTrimmed-mean(costTrimmed))>costSTD)
costTrimmed(abs(costTrimmed-mean(costTrimmed))>costSTD) = [];
plot(1:length(costTrimmed),costTrimmed)
title('Cost Over Time')

ylabel('Cost')


