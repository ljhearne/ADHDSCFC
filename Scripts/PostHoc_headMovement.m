% Hsiang-Yuan Lin
% Attachments
% Sun, Jun 30, 3:37 AM (1 day ago)
% to me, Luca
% 
% Dear Luke,
% 
% In the .mat, 
% ADHD_M(:,5)=mean FD
% ADHD_M(:,6)=DSI signal dropout counts (motion estimates)
% ADHD_M(:,7)=group assignment, 1=CTRL, 2=ADHD
% 
% score(:,1)= factor scores of PCA component 1 based on ADHD_M(:,[1:4])
% 
% Cheers, Hsiang-Yuan

clearvars
close all

load('ADHD_PCA_Corr.mat','ADHD_M','score')

fd = ADHD_M(:,5);
dsi_motion = ADHD_M(:,6);
group = ADHD_M(:,7);
symptoms = score(:,1);

%% print stats to screen
group_label = {'ADHD','CTRL'};
for i = 1:2
    [r,p] = corr(fd(group==i),symptoms(group==i));
    disp(['Measure: framewise displacement, Group: ',group_label{i},', r=',num2str(r),', p=',num2str(p)])
    
    [r,p] = corr(dsi_motion(group==i),symptoms(group==i));
    disp(['Measure: dsi motion, Group: ',group_label{i},', r=',num2str(r),', p=',num2str(p)])
    
end

%% plots
addpath(genpath('Functions'));
[cb] = cbrewer('qual','Set3',12,'pchip');
cl(1,:) = [0.5 0.5 0.5];
cl(2,:) = cb(1,:);

figure('Color','w','Position',[1200 425 400 150]); hold on
subplot(1,2,1)
for i = 1:2
    scatter(fd(group==i),symptoms(group==i),...
        'MarkerEdgeColor',cl(i,:),...
        'MarkerFaceColor',cl(i,:),...
        'MarkerEdgeAlpha', 0.7,...
        'MarkerFaceAlpha', 0.7,...
        'SizeData',20,...
        'LineWidth',1); hold on;
    
    h = lsline;     
    set(h(1),'LineWidth',2,'Color',cl(i,:)*0.6);
    hold on
end

set(gca,'FontName', 'Helvetica','FontSize', 10,'box','off');
ylabel('ADHD Symptoms ');
xlabel('Framewise displacement');

subplot(1,2,2)
for i = 1:2
    scatter(dsi_motion(group==i),symptoms(group==i),...
        'MarkerEdgeColor',cl(i,:),...
        'MarkerFaceColor',cl(i,:),...
        'MarkerEdgeAlpha', 0.7,...
        'MarkerFaceAlpha', 0.7,...
        'SizeData',20,...
        'LineWidth',1); hold on;
    
    h = lsline;     
    set(h(1),'LineWidth',2,'Color',cl(i,:)*0.6);
    hold on
end

set(gca,'FontName', 'Helvetica','FontSize', 10,'box','off');
ylabel('ADHD Symptoms ');
xlabel('DSI motion');

%save figure
saveas(gcf,'/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/Figures/head_motion_plot.jpeg');