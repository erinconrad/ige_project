function log_rank_erin(x1,x2,results_folder,feature,duration)

%{
This function is called by the main function ige_stats and generates the
Kaplan Meier curve for the survival analysis
%}

all_x = [x1;x2];
J = max(all_x(:,1)); % max time
x_cell = {x1,x2};
starting_num = [size(x1,1),size(x2,1)];



surv_prob = ones(J,2);
int_surv_prob = ones(J,2);


N = ones(J,2).*starting_num;
E = ones(J,2);
O = zeros(J,2);
V = ones(J,2);

for j = 2:J
    for i = 1:2
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


%% Get censoring times
cens = cell(2,1);
for i = 1:2
    x = x_cell{i};
    cens{i} = x(x(:,2) == 1,1);
end

%% Kaplan Meier plot
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
l1=plot([1:J]/60,surv_prob(:,1),'k','linewidth',2); % percentage of surviving things over time for each group
hold on
l2=plot([1:J]/60,surv_prob(:,2),'k--','linewidth',2);
c = plot(cens{1}/60,surv_prob(cens{1},1),'+','color',[0, 0.4470, 0.7410],'MarkerSize',20,'linewidth',2);
plot(cens{2}/60,surv_prob(cens{2},2),'+','color',[0, 0.4470, 0.7410],'MarkerSize',20,'linewidth',2)
legend([l1,l2,c],{'Drug resistant','Drug responsive','End of EEG'},'location','southeast')
ylim([0 1])
xlabel('Time (hours)')
ylabel(sprintf('Probability of not having %s',feature))
set(gca,'fontsize',20)
annotation('textbox',[0.03 0.45 0.1 0.1],'String','B',...
    'linestyle','none','fontsize',35);
print(gcf,[results_folder,'surv_plot_',feature],'-depsc')






end