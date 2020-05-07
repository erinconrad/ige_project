function log_rank_erin(x_cell,results_folder,legend_labels,feature,duration)

%{
This function is called by the main function ige_stats and generates the
Kaplan Meier curve for the survival analysis
%}


%% Get the maximum time
% Max time per thing
max_per_thing = cellfun(@(x) max(x(:,1)),x_cell,'UniformOutput',false);
max_per_thing = cell2mat(max_per_thing);
J = max(max_per_thing); % max across all things

starting_num = cellfun(@(x) size(x,1),x_cell,'UniformOutput',false);
starting_num = cell2mat(starting_num);


surv_prob = ones(J,length(x_cell));
int_surv_prob = ones(J,length(x_cell));


N = ones(J,length(x_cell)).*starting_num;
E = ones(J,length(x_cell));
O = zeros(J,length(x_cell));
V = ones(J,length(x_cell));

for j = 2:J
    for i = 1:length(x_cell)
        x = x_cell{i};
 
        % Number who have "died" before this time
        n_died = sum(x(:,1) < j & x(:,2) == 0);
        
        % Number censored before this time
        n_cens = sum(x(:,1) < j & x(:,2) == 1);
        
        % Get number at risk (starting number minus those censored or dead)
        n_at_risk = starting_num(i) - n_died - n_cens;
        
        % Number who died this interval
        n_died_this_interval = sum(x(:,1) == j & x(:,2) == 0);

        % Number who survived this interval
        n_surv = n_at_risk - n_died_this_interval;
        
        % Interval survival probability
        int_surv_prob(j,i) = n_surv/n_at_risk;
        
        % Surv_prob
        surv_prob(j,i) = prod(int_surv_prob(1:j,i));
                
        % Stuff for test
        N(j,i) = n_at_risk;
        O(j,i) = n_died_this_interval;
        
    end
end

%{
%% Stuff for my log-rank test (doesn't work; I use R for this)

% Loop again to get E and V
for j = 1:J
    for i = 1:2

        E(j,i) = N(j,i)*(O(j,1)+O(j,2))/(N(j,1)+N(j,2));
        V(j,i) = E(j,i)*(N(j,1)+N(j,2)-(O(j,1)+O(j,2)))/(N(j,1)+N(j,2))*...
            (N(j,1)+N(j,2)-N(j,i))/(N(j,1)+N(j,2)-1);
    end
end

% Remove nans
E(isnan(E)) = 0;
V(isnan(V)) = 0;

Z_num = 0;

% define i to be 1
i = 1;

% Loop to get Z
for j = 1:J
    Z_num = Z_num + (O(j,1) - E(j,i));
end

Z_sq_denom = 0;
for j = 1:J
    Z_sq_denom = Z_sq_denom + V(j,i);
end

Z_denom = sqrt(Z_sq_denom);

Z = Z_num/Z_denom;
p=2*(1-0.5*erfc(-Z/realsqrt(2)));

%}

%% Get censoring times
cens = cell(length(x_cell),1);
for i = 1:length(x_cell)
    x = x_cell{i};
    cens{i} = x(x(:,2) == 1,1);
end

%% Kaplan Meier plot
if length(x_cell) == 2
    figure
    set(gcf,'position',[215 700 900 600]);
    pl1 = subplot(2,1,1);
    pl1_pos = get(pl1,'position');
    %set(pl1,'position',[pl1_pos(1) pl1_pos(2) pl1_pos(3)*2/3 pl1_pos(4)])
    histogram(duration/60,30,'facecolor',[0, 0.4470, 0.7410],'linewidth',2)
    xlabel('Total EEG duration (hours)')
    ylabel('Number of patients')
    set(gca,'fontsize',20)
    annotation('textbox',[0.015 0.88 0.1 0.1],'String','A',...
        'linestyle','none','fontsize',35);

    subplot(2,1,2)
    
    l1=plot([1:J]/60,surv_prob(:,1),'k','linewidth',2); % percentage of surviving things over time for each group
    hold on
    l2=plot([1:J]/60,surv_prob(:,2),'k--','linewidth',2);
    c = plot(cens{1}/60,surv_prob(cens{1},1),'+','color',[0, 0.4470, 0.7410],'MarkerSize',20,'linewidth',2);
    plot(cens{2}/60,surv_prob(cens{2},2),'+','color',[0, 0.4470, 0.7410],'MarkerSize',20,'linewidth',2)
    legend([l1,l2,c],{'Drug resistant','Drug responsive','End of EEG'},'location','southwest')
    ylim([0 1])
    xlabel('Time (hours)')
    ylabel(sprintf('Probability of not having %s',feature))
    set(gca,'fontsize',20)
    annotation('textbox',[0.015 0.45 0.1 0.1],'String','B',...
        'linestyle','none','fontsize',35);
    print(gcf,[results_folder,'Figure2'],'-depsc')
else
    
    colors = [0, 0.4470, 0.7410;...
        0.8500, 0.3250, 0.0980;...
        0.9290, 0.6940, 0.1250;...
        0.4940, 0.1840, 0.5560];
    
    figure
    set(gcf,'position',[215 700 900 600]);
    pl1 = subplot(2,1,1);
    pl1_pos = get(pl1,'position');
    %set(pl1,'position',[pl1_pos(1) pl1_pos(2) pl1_pos(3)*2/3 pl1_pos(4)])
    histogram(duration/60,30,'facecolor',[0, 0.4470, 0.7410],'linewidth',2)
    xlabel('Total EEG duration (hours)')
    ylabel('Number of patients')
    set(gca,'fontsize',20)
    annotation('textbox',[0.03 0.88 0.1 0.1],'String','A',...
        'linestyle','none','fontsize',35);

    subplot(2,1,2)
    l = zeros(length(x_cell),1);
    for i = 1:length(x_cell)
        l(i) = plot([1:J]/60,surv_prob(:,i),'color',colors(i,:),'linewidth',2);
        hold on
        plot(cens{i}/60,surv_prob(cens{i},i),'+','color',colors(i,:),'MarkerSize',20,'linewidth',2);
    end
    %xl = get(gca,'xlim');
    c = plot(inf,inf,'k+','MarkerSize',20,'linewidth',2);
    legend_stuff = [l;c];
    legend(legend_stuff,legend_labels,'location','southeast')
    ylim([0 1])
    xlabel('Time (hours)')
    ylabel(sprintf('Probability of not having %s',feature))
    set(gca,'fontsize',20)
    annotation('textbox',[0.03 0.45 0.1 0.1],'String','B',...
        'linestyle','none','fontsize',35);
    print(gcf,[results_folder,'surv_plot_',feature],'-depsc')
        
end






end