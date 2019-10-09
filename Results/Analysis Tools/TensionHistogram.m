clc
clear all

%samples will be set up as column vectors
name = '../2018-12-24_RollingXDirectionPairedCables_slowMotors';
file_list = dir([name,'/*.mat']);
file_list = {file_list.name};
tensionList = [];
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
    
    tensionList = [tensionList,tensions];
end

%% plots
figure(1)
tensions = sum(tensionList);
histogram(tensions,50,'normalization','probability')
hold on
initialtension = tensions(1);
line([initialtension initialtension],[0 0.2],'linewidth',4,'LineStyle','--','Color','y')
xlabel('Instantaneous Total Cable Tension  [N]')
ylabel('Relative Probability')
avgTension = mean(tensions);
line([avgTension avgTension],[0 0.2],'linewidth',4,'LineStyle','--','Color','g')
set(gca,'fontsize',20)
set(gca,'fontname','times new roman')
