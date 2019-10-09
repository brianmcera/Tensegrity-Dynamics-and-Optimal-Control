clc
clear all

%samples will be set up as column vectors
name = '../2019-03-22_PairedCable_fastMotors_xDirection_perturbedCables';
file_list = dir([name,'/*.mat']);
file_list = {file_list.name};
XspeedList = [];
YspeedList = [];
ZspeedList = [];
for file_num = 1:numel(file_list)
    filename = [name,'/',file_list{file_num}];
    disp(['Loading ',filename])
    try
        load(filename)
    catch
        disp(['Error loading file: ',filename])
        continue %if there's errors loading file, skip file
    end
    t_time = find((X_record.p(1,:)),1,'last');
    XspeedList = [XspeedList,mean(X_record.pDOT(1:3:end,1:t_time))];
    YspeedList = [YspeedList,mean(X_record.pDOT(2:3:end,1:t_time))];
    ZspeedList = [ZspeedList,mean(X_record.pDOT(3:3:end,1:t_time))];
end
totalSpeedList = mean([XspeedList;YspeedList;ZspeedList]);

%% sort CoM Speed
averageSpeedList = zeros(1,numel(file_list))
for file_num = 1:numel(file_list)
    filename = [name,'/',file_list{file_num}];
    disp(['Loading ',filename])
    try
        load(filename)
    catch
        disp(['Error loading file: ',filename])
        continue %if there's errors loading file, skip file
    end
    t_start = 1;
    dt = simulationParameters_record.timestep;
    t_time = find((X_record.p(1,:)),1,'last');
    COM_x = mean(X_record.p(1:3:end,:));
    COM_y = mean(X_record.p(2:3:end,:));
    averageSpeedList(file_num) = norm([COM_x(t_time),COM_y(t_time)]-...
        [COM_x(t_start),COM_y(t_start)])/(dt*(t_time-t_start));
end

%% plots
figure(1)
histogram(XspeedList,50,'normalization','probability')
hold on
line([.018 .018],[0 0.1],'linewidth',4,'LineStyle','--','Color','r')
xlabel('Instantaneous Speed of Center of Mass [m/s]')
ylabel('Relative Probability')
meanXSpeed = mean(XspeedList);
line([meanXSpeed meanXSpeed],[0 0.1],'linewidth',4,'LineStyle','--','Color','g')
set(gca,'fontsize',20)
set(gca,'fontname','times new roman')
