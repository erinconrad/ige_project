function parameter = ige_stats

%% Parameters

% minimum number of elements in any cell of contingency table to do chi2
min_num_fisher = 20; % This is very conservative

% Establish relevant parameters
parameter(1).name = 'drug_resistant';
parameter(2).name = 'eeg_duration';
parameter(3).name = 'sex';
parameter(4).name = 'sleep';
parameter(5).name = 'vpa';
parameter(6).name = 'exclude_eeg___0';
parameter(7).name = 'gsw___1';
parameter(8).name = 'gsw___2';
parameter(9).name = 'psw___1';
parameter(10).name = 'psw___2';
parameter(11).name = 'pst___1';
parameter(12).name = 'pst___2';
parameter(13).name = 'gpfa___1';
parameter(14).name = 'gpfa___2';
parameter(15).name = 'glvfa___1';
parameter(16).name = 'glvfa___2';
parameter(17).name = 'slow___1';
parameter(18).name = 'slow___2';
parameter(19).name = 'foc_dis___1';
parameter(20).name = 'foc_dis___2';
parameter(21).name = 'foc_slow___1';
parameter(22).name = 'foc_slow___2';
parameter(23).name = 'exclude_clinical___0';
parameter(24).name = 'duration_minutes'; % this is not in the original table; it will be a new column
parameter(25).name = 'gsw';
parameter(26).name = 'psw';
parameter(27).name = 'pst';
parameter(28).name = 'gpfa';
parameter(29).name = 'glvfa';
parameter(30).name = 'slow';
parameter(31).name = 'foc_dis';
parameter(32).name = 'foc_slow';
parameter(33).name = 'age_at_eeg';

% tell it which parameters are the awake version of the one below it
eeg_parameters = [7:22,25:32];
awake_parameters = [7:2:21];
summary_parameters = [25:32];

%% File path
csv_path = "/Users/erinconrad/Desktop/residency stuff/R25/ige project/data/IGEDatabase-DeidentifiedData_DATA_2020-02-20_2017.csv";
r_file_path = "/Users/erinconrad/Desktop/residency stuff/R25/ige project/data/data_for_r.csv";
results_folder = '/Users/erinconrad/Desktop/residency stuff/R25/ige project/results/';

%% Load csv file
data = readtable(csv_path);

%% Get column numbers
new_col = 1;
for i = 1:length(parameter)
    parameter(i).column = find(strcmp(data.Properties.VariableNames,parameter(i).name));

    if isempty(parameter(i).column) == 1
        if i >= 24
            parameter(i).column = size(data,2)+new_col;
            new_col = new_col + 1;
        else
            error('what\n');
        end
    end
end
feature_cols = parameter(7).column:parameter(22).column;
age_eeg_col = parameter(33).column;


%% Prep new table
new_table = data;

% Find eeg rows in table
eeg_row = zeros(size(data,1),1);
curr_num = 0;
for i = 1:size(data,1)
    if curr_num == data.record_id(i) % if it's the same record id as the last one, it's an eeg row
        eeg_row(i) = 1;
    end
    curr_num = data.record_id(i);
end

% Prep array for duration in minutes
duration_minutes = zeros(size(eeg_row));

% Prep arrays for new eeg summary features (concatenating sleep and
% wakefulness)
for i = 1:length(summary_parameters)
    summ(i).array = zeros(size(eeg_row));
end

% Loop through rows and add eeg data to main clinical row
for i = 1:size(data,1)-1
    
    % continue if it's an eeg row
    if eeg_row(i) == 1, continue; end
    
    % Find the rows corresponding to the eeg data for this clinical row
    curr_eeg_rows = [];
    k = 1; % the possible eeg row we're checking
    while 1
        if i+k >size(data,1)
            break % break if it's end of table
        end
        if eeg_row(i+k) == 1 % add it if it's an eeg row
            curr_eeg_rows = [curr_eeg_rows;i+k];
        else
            break % once you hit the first non-eeg row, you're done
        end
        k = k + 1; % move to the next row
    end
    
    % Get if there is any stage 2 sleep in any eeg
    sleep = table2array(data(curr_eeg_rows,parameter(4).column));
    if any(sleep == 1) % Note that 1 means stage 2 sleep included
        new_table.sleep(i) = 1;
    else
        new_table.sleep(i) = 0;
    end
    
    % Get total duration of eeg
    duration1 = table2array(data(curr_eeg_rows,parameter(2).column));
    total_dur = 0;
    for j = 1:length(duration1)
        dur_text = duration1{j};
        
        % Fix for mistakes in text
        dur_text = strrep(dur_text,'.',':');
        
        if length(dur_text) == 4 || length(dur_text) == 5
            dur_num = duration(dur_text,'InputFormat','hh:mm');
        elseif length(dur_text) == 8
            dur_num = duration(dur_text,'InputFormat','hh:mm:ss');
        elseif isempty(dur_text) == 1
            dur_num = 0;
        elseif strcmp(dur_text,'41 minutes') == 1 % correction for mistake in text
            dur_num = duration('00:41','InputFormat','hh:mm');
        else
            fprintf('Warning, duration text for row %d says %s\n',i,dur_text);
            dur_num = 0;
        end
        total_dur = dur_num + total_dur;

    end
    
    % convert duration to minutes to get a single number
    if total_dur == 0
        duration_minutes(i) = 0;
    else
        duration_minutes(i) = minutes(total_dur);
    end
    
    % summarize eeg features
    for j = feature_cols
        curr_col = table2array(data(curr_eeg_rows,j));
        
        % If any are checked, make the clinical column be 1
        if any(curr_col) == 1
            new_table{i,j} = 1;
        else
            new_table{i,j} = 0;
        end
    end

    % Get age at first eeg
    curr_col = table2array(data(curr_eeg_rows,age_eeg_col));
    new_table{i,age_eeg_col} = min(curr_col);
    
    % now concatenate sleep and wake for new columns
    count = 0;
    for j = awake_parameters
        count = count + 1; % which of the awake parameters it is
        
        % get the columns corresponding to this AND the sleep version
        curr_cols = table2array(data(curr_eeg_rows,...
            parameter(j).column:parameter(j+1).column));
        
        % for testing
        if 0
        fprintf('\n\n');
        data.record_id(i)
        parameter(j).name
        curr_cols
        any(curr_cols(:))
        pause
        end
        
        % If any are checked in awake OR sleep, make the new clinical
        % column be 1
        if any(curr_cols(:)) == 1
            summ(count).array(i) = 1;
        else
            summ(count).array(i) = 0;
        end
        
    end
    
    % Get if we are excluding based on the eeg
    eeg_exclude = table2array(data(curr_eeg_rows,parameter(6).column));
    if any(eeg_exclude == 1)
        new_table.exclude_eeg___0(i) = 1;
    else
        new_table.exclude_eeg___0(i) = 0;
    end
    
end


% Add duration info
new_table = addvars(new_table,duration_minutes);

% Add summary info
for i = 1:length(summ)
    new_table = addvars(new_table,summ(i).array);
end

% Rename summary variable names
for i = summary_parameters
    new_table.Properties.VariableNames{parameter(i).column} = parameter(i).name;
end


% Remove eeg rows
new_table(logical(eeg_row),:) = [];


%% Identify and remove exclusion rows
exclude = any([new_table.exclude_clinical___0==1,...
    new_table.exclude_eeg___0 == 1],2);

% Remove exclusion rows
new_table(exclude,:) = [];

fprintf('Excluded %d EEGs per criteria.\n',sum(exclude));

%% Identify and remove zero duration eegs (ideally shouldn't have these)
zero_duration = find(new_table.duration_minutes == 0);
for i = 1:length(zero_duration)
    fprintf('\nWarning, record ID %d with zero duration. Removing...\n',...
        new_table.record_id(zero_duration(i)));
end
new_table(zero_duration,:) = [];

%% Redefine drug resistance to be 1 or 0
% In the Redcap table, drug_resistant = 2 means responsive, 1 means
% resistant, so now 1 means responsive, 0 means resistant
new_table.drug_resistant = new_table.drug_resistant - 1;

% So flip it so now 1 means resistant, 0 means responsive
new_table.drug_resistant = ~new_table.drug_resistant;

%% See if any nans for drug resistant (there shouldn't be)
dr_nan = find(isnan(new_table.drug_resistant));
if isempty(dr_nan) == 0
    for i = 1:length(dr_nan)
        fprintf('\nDrug resistance is nan for record id %d.\n',new_table.record_id(dr_nan(i)));
    end
end

%% Print some random rows to test that I didn't mess up
if 0
    for i = 1:10
        r = randi(size(new_table,1));
        new_table(r,:)
    end
end

%{
********************************
Now do some statistics
********************************
%}

%% Wilcoxon rank sums for duration and age at first eeg
for pdur = [24 33]
    name = parameter(pdur).name;
    [pval,~,stats] = ranksum(new_table.(name)(new_table.drug_resistant == 1),...
        new_table.(name)(new_table.drug_resistant == 0));
    [pval_alt,~,U] = ranksum_output_u(new_table.(name)(new_table.drug_resistant == 1),...
        new_table.(name)(new_table.drug_resistant == 0));
    
    if pval ~= pval_alt, error('what\n'); end
    
    parameter(pdur).stats.U = U;
    parameter(pdur).stats.p = pval;
    parameter(pdur).stats.other = stats;
    parameter(pdur).stats.avg_dur = [nanmean(new_table.(name)(new_table.drug_resistant == 1)),...
        nanmean(new_table.(name)(new_table.drug_resistant == 0))];
    fprintf(['\nAverage %s is %1.1f m for drug resistant patients and \n'...
        '%1.1f m for drug responsive patients (Wilcoxon rank sum: p = %1.3f).\n'],...
        name,parameter(pdur).stats.avg_dur(1),parameter(pdur).stats.avg_dur(2),parameter(pdur).stats.p);
    
end

% Also loop through and do WRSs for other parameters
duration_feature_table = cell2table(cell(0,5),'VariableNames',{'Feature',...
    'AverageDurationWithFeature','AverageDurationWithoutFeature','U','p'});
for p = 25:32
    % get appropriate column
    name = parameter(p).name;
    
    % Get number who have parameter
    par = new_table.(name);
    
    if sum(par) == 0, continue; end
    
    [pval,~,stats] = ranksum(new_table.duration_minutes(par == 1),...
    new_table.duration_minutes(par == 0));
    [alt_pval,~,U] = ranksum_output_u(new_table.duration_minutes(par == 1),...
    new_table.duration_minutes(par == 0));

    if isequal(pval,alt_pval) == 0, error('what\n'); end

    duration_feature_table = [duration_feature_table;...
        cell2table({name,nanmean(new_table.duration_minutes(par == 1)),...
        nanmean(new_table.duration_minutes(par == 0)),U,pval},'VariableNames',{'Feature',...
    'AverageDurationWithFeature','AverageDurationWithoutFeature','U','p'})];

    if 0
        fprintf(['\nAverage duration is %1.1f m for patients with %s and \n'...
        '%1.1f m for patients without %s (Wilcoxon rank sum: p = %1.3f).\n'],...
        nanmean(new_table.duration_minutes(par == 1)),name,...
        nanmean(new_table.duration_minutes(par == 0)),name,pval);
    end
end

% Reformat table
duration_feature_table.Feature{1} = 'GSW';
duration_feature_table.Feature{2} = 'PSW';
duration_feature_table.Feature{3} = 'PST';
duration_feature_table.Feature{4} = 'GPFA';
duration_feature_table.Feature{5} = 'GLVFA';
duration_feature_table.Feature{6} = 'Focal discharges';
duration_feature_table.Feature{7} = 'Focal slowing';
duration_feature_table.AverageDurationWithFeature = ...
    arrayfun(@(x)sprintf('%1.1f',x),duration_feature_table.AverageDurationWithFeature,'UniformOutput',false);
duration_feature_table.AverageDurationWithoutFeature = ...
    arrayfun(@(x)sprintf('%1.1f',x),duration_feature_table.AverageDurationWithoutFeature,'UniformOutput',false);
duration_feature_table.p = arrayfun(@(x)sprintf('%1.3f',x),duration_feature_table.p,'UniformOutput',false);

%duration_feature_table.Properties.VariableNames{'AverageDurationWithFeature'} = 'Duration with feature';
%duration_feature_table.Properties.VariableNames{'AverageDurationWithoutFeature'} = 'Duration without feature';
%duration_feature_table.Properties.VariableNames{'p'} = 'p value';


if 0
%% Figure for percentage of patients with feature at different durations
% X axis is duration of eeg
% Y axis is number of patients whose eeg is that duration or lower who have
% the feature of interest
fig1 = figure;
set(fig1,'position',[100 416 600 384])

fig2 = figure;
set(fig2,'position',[100 416 600 384])
legend_names = {};

features_to_plot = 25:29;

dur_short = zeros(length(features_to_plot),1);
count = 0;
for p = features_to_plot
    count = count + 1;
    
    % get appropriate column
    name = parameter(p).name;
    legend_names = [legend_names,name];
    
    % Get number who have parameter
    par = new_table.(name);
    
    % Get duration
    dur = new_table.duration_minutes;
    
    % Define duration step size
    dur_steps = 0:20:max(dur);
    n_par_dur = zeros(1,length(dur_steps));
    perc_par_dur = zeros(1,length(dur_steps));
    
    % Loop through duration steps
    for d = 1:length(dur_steps)
        
        % Get the number of patients who have parameter and have duration
        % less than or equal to specified duration
        n_par_dur(d) = sum(par(dur<=dur_steps(d)));
        
        % Get the percentage of patients who have parameter within that
        % specified duration (divide by the total number with that
        % parameter)
        perc_par_dur(d) = sum(par(dur<=dur_steps(d)))/sum(par)*100;
    end
    figure(fig1)
    plot(dur_steps/60,perc_par_dur,'linewidth',2)
    hold on
    
    % Compare less than 23 hour EEGs
    dur_short(count) = sum(par(dur<=12*60))/sum(par)*100;
    
end

figure(fig1)
xlabel('EEG duration (hours)')
ylabel('Percentage of patients with feature');
legend(legend_names,'location','southeast')
set(gca,'fontsize',20)
print(fig1,[results_folder,'Figure1'],'-depsc')
close(fig1)

figure(fig2)
bar(dur_short)
xticklabels(legend_names)
set(gca,'fontsize',20)
xlabel('EEG feature')
ylabel({'Percentage of patients with feature';'who have it in <24 hours'});
end


%% Categoricla Univariate stats
% Loop through categorical parameters of interest
for p = [3:5,eeg_parameters]
    
    % Note that male == 1, female == 0
    % Treid VPA == 1, not tried VPA == 0
    
    % get appropriate column
    name = parameter(p).name;
    
    % Get number who have parameter
    par = new_table.(name);
    
    % do chi squared
    [tbl,chi2,pval,labels] = crosstab(new_table.drug_resistant,...
        par);
    
    % Don't do the analysis if all par is 0
    if length(unique(par)) == 1
        parameter(p).stats.tbl = [tbl,[0;0]];
        parameter(p).stats.other = nan;
        continue
    end
    
    % if the value of any cell in the contingency table is <=min_num_fisher, do fisher
    % exact test instead
    if any(tbl(:)) <= min_num_fisher
        [~,pval,stats] = fishertest(tbl);
        parameter(p).stats.labels = labels;
        parameter(p).stats.tbl = tbl;
        parameter(p).stats.p = pval;
        parameter(p).stats.OddsRatio = stats.OddsRatio;
        parameter(p).stats.ConfidenceInterval = stats.ConfidenceInterval;
    else
        fprintf('\nNote that I am doing chi2 test for %s.\n',parameter(p).name);
        parameter(p).stats.tbl = tbl;
        parameter(p).stats.chi2 = chi2;
        parameter(p).stats.p = pval;
        parameter(p).stats.labels = labels;
    end
    
    
    
    
end

%% Table 1
% This will be long and ugly

% Total number
clinical_table = cell2table({'TotalNumber',sprintf('%d (%1.1f%%)',...
    sum(new_table.drug_resistant == 0),...
    sum(new_table.drug_resistant == 0)/(sum(new_table.drug_resistant == 0)+sum(new_table.drug_resistant == 1))*100),...
    sprintf('%d (%1.1f%%)',...
    sum(new_table.drug_resistant == 1),...
    sum(new_table.drug_resistant == 1)/(sum(new_table.drug_resistant == 0)+sum(new_table.drug_resistant == 1))*100),'',''},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'});

% Sex
clinical_table = [clinical_table;cell2table({'Sex','','','',''},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];
p = 3;

clinical_table = [clinical_table;cell2table({'Men',...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(1,2),...
    parameter(p).stats.tbl(1,2)/(parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(2,2))*100),...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(2,2),...
    parameter(p).stats.tbl(2,2)/(parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(2,2))*100),...
    sprintf('%1.2f',parameter(p).stats.OddsRatio),...
    sprintf('%1.3f',parameter(p).stats.p)},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];

clinical_table = [clinical_table;cell2table({'Women',...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(1,1),...
    parameter(p).stats.tbl(1,1)/(parameter(p).stats.tbl(1,1)+parameter(p).stats.tbl(2,1))*100),...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(2,1),...
    parameter(p).stats.tbl(2,1)/(parameter(p).stats.tbl(1,1)+parameter(p).stats.tbl(2,1))*100),...
    '',''},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];


% Age at first eeg
p = 33;
clinical_table = [clinical_table;cell2table({parameter(p).name,...
    sprintf('%1.1f',parameter(p).stats.avg_dur(2)),sprintf('%1.1f',parameter(p).stats.avg_dur(1))...
    sprintf('%1.2f',parameter(p).stats.U),sprintf('%1.3f',parameter(p).stats.p)},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];

% Tried VPA
clinical_table = [clinical_table;cell2table({'TriedVPA','','','',''},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];
p = 5;

clinical_table = [clinical_table;cell2table({'Yes',...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(1,2),...
    parameter(p).stats.tbl(1,2)/(parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(2,2))*100),...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(2,2),...
    parameter(p).stats.tbl(2,2)/(parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(2,2))*100),...
    sprintf('%1.2f',parameter(p).stats.OddsRatio),sprintf('%1.3f',parameter(p).stats.p)},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];

clinical_table = [clinical_table;cell2table({'No',...
    sprintf('%d (%1.1f%%)',...    
    parameter(p).stats.tbl(1,1),...
    parameter(p).stats.tbl(1,1)/(parameter(p).stats.tbl(1,1)+parameter(p).stats.tbl(2,1))*100),...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(2,1),...
    parameter(p).stats.tbl(2,1)/(parameter(p).stats.tbl(1,1)+parameter(p).stats.tbl(2,1))*100),...
    '',''},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];

% EEG duration
p = 24;
clinical_table = [clinical_table;cell2table({parameter(p).name,...
    sprintf('%1.1f',parameter(p).stats.avg_dur(2)),...
    sprintf('%1.1f',parameter(p).stats.avg_dur(1)),...
    sprintf('%1.2f',parameter(p).stats.U),sprintf('%1.3f',parameter(p).stats.p)},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];

% Captured sleep
clinical_table = [clinical_table;cell2table({'CapturedSleep','','','',''},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];
p = 4;

clinical_table = [clinical_table;cell2table({'Yes',...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(1,2),...
    parameter(p).stats.tbl(1,2)/(parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(2,2))*100),...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(2,2),...
    parameter(p).stats.tbl(2,2)/(parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(2,2))*100),...
    sprintf('%1.2f',parameter(p).stats.OddsRatio),...
    sprintf('%1.3f',parameter(p).stats.p)},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];

clinical_table = [clinical_table;cell2table({'No',...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(1,1),...
    parameter(p).stats.tbl(1,1)/(parameter(p).stats.tbl(1,1)+parameter(p).stats.tbl(2,1))*100),...
    sprintf('%d (%1.1f%%)',...
    parameter(p).stats.tbl(2,1),...
    parameter(p).stats.tbl(2,1)/(parameter(p).stats.tbl(2,1)+parameter(p).stats.tbl(1,1))*100),...
    '',''},...
    'VariableNames',{'Parameter','Responsive','Resistant','Statistic','P'})];

% Format columns nicely
clinical_table.Parameter{1} = 'Total number';
clinical_table.Parameter{5} = 'Age at first EEG';
clinical_table.Parameter{6} = 'Tried VPA?';
clinical_table.Parameter{9} = 'Duration (minutes)';
clinical_table.Parameter{10} = 'Captured sleep?';
%clinical_table.Responsive = cellfun(@(x)sprintf('%1.1f',x),clinical_table.Responsive,'UniformOutput',false);
%clinical_table.Resistant = cellfun(@(x)sprintf('%1.1f',x),clinical_table.Resistant,'UniformOutput',false);
%clinical_table.Statistic = cellfun(@(x)sprintf('%1.1f',x),clinical_table.Statistic,'UniformOutput',false);
%clinical_table.P = cellfun(@(x)sprintf('%1.3f',x),clinical_table.P,'UniformOutput',false);

fprintf('\n\nTable 1:\n');
clinical_table
fprintf('\n');
writetable(clinical_table,[results_folder,'Table1.csv']);


fprintf('\n\nTable2:\n');
duration_feature_table
fprintf('\n');
writetable(duration_feature_table,[results_folder,'Table2.csv']);

%% Make a summary table of univariate statistics for eeg features
all_names = {};
all_p_value = [];
all_stat = [];
all_n_resistant = [];
all_n_responsive = [];

all_perc_resistant = [];
all_perc_responsive = [];

for p = 25:32
    if strcmp(parameter(p).name,'duration_minutes') == 1, continue; end
    if strcmp(parameter(p).name,'sex') == 1, continue; end
    
    if isfield(parameter(p).stats,'p') == 0
        continue
    end
    
    all_names = [all_names;parameter(p).name];
    all_p_value = [all_p_value;parameter(p).stats.p];
    all_stat = [all_stat;parameter(p).stats.OddsRatio];
    %{
    if isfield(parameter(p).stats,'OddsRatio')
        all_stat = [all_stat;parameter(p).stats.OddsRatio];
    elseif isfield(parameter(p).stats,'chi2')
        all_stat = [all_stat;parameter(p).stats.chi2];
    end
    %}
    
    all_n_resistant = [all_n_resistant;parameter(p).stats.tbl(2,2)];
    all_n_responsive = [all_n_responsive;parameter(p).stats.tbl(1,2)];
    
    all_perc_resistant = [all_perc_resistant;parameter(p).stats.tbl(2,2)/...
        (parameter(p).stats.tbl(2,2)+parameter(p).stats.tbl(2,1))];
    all_perc_responsive = [all_perc_responsive;parameter(p).stats.tbl(1,2)/...
        (parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(1,1))];
    
end

str_resistant = cell(length(all_n_resistant),1);
str_responsive = cell(length(all_n_resistant),1);

for i = 1:length(all_n_resistant)
    str_resistant{i} = sprintf('%d (%1.1f%%)',all_n_resistant(i),all_perc_resistant(i)*100);
    str_responsive{i} = sprintf('%d (%1.1f%%)',all_n_responsive(i),all_perc_responsive(i)*100);
end

eeg_feature_table = table(all_names,str_responsive,str_resistant,all_stat,all_p_value,...
    'VariableNames',{'Feature','Responsive','Resistant',...
    'OddsRatio','PValue'});

% Format table
eeg_feature_table.OddsRatio = arrayfun(@(x)sprintf('%1.2f',x),eeg_feature_table.OddsRatio,'UniformOutput',false);
eeg_feature_table.PValue = arrayfun(@(x)sprintf('%1.3f',x),eeg_feature_table.PValue,'UniformOutput',false);
eeg_feature_table.Feature{1} = 'GSW';
eeg_feature_table.Feature{2} = 'PSW';
eeg_feature_table.Feature{3} = 'PST';
eeg_feature_table.Feature{4} = 'GPFA';
eeg_feature_table.Feature{5} = 'GLVFA';
eeg_feature_table.Feature{6} = 'Focal discharges';
eeg_feature_table.Feature{7} = 'Focal slowing';
%eeg_feature_table.Properties.VariableNames{'OddsRatio'} = 'Odds ratio';
%eeg_feature_table.Properties.VariableNames{'PValue'} = 'p value';

fprintf('\n\nTable3:\n');
eeg_feature_table
fprintf('\n');
writetable(eeg_feature_table,[results_folder,'Table3.csv']);


%% Supplemental table for sleep vs wake features
all_names = {};
all_p_value = [];
all_stat = [];
all_n_resistant = [];
all_n_responsive = [];

all_perc_resistant = [];
all_perc_responsive = [];

for p = [7:16,19:22] %wake and sleep eeg features
    if strcmp(parameter(p).name,'duration_minutes') == 1, continue; end
    if strcmp(parameter(p).name,'sex') == 1, continue; end
    
    
    
    all_names = [all_names;parameter(p).name];
    if isfield(parameter(p).stats,'p') == 0
        all_p_value = [all_p_value;nan];
        all_stat = [all_stat;nan];
    else
        all_p_value = [all_p_value;parameter(p).stats.p];
        all_stat = [all_stat;parameter(p).stats.OddsRatio];
    end
    
    %{
    if isfield(parameter(p).stats,'OddsRatio')
        all_stat = [all_stat;parameter(p).stats.OddsRatio];
    elseif isfield(parameter(p).stats,'chi2')
        all_stat = [all_stat;parameter(p).stats.chi2];
    end
    %}

    all_n_resistant = [all_n_resistant;parameter(p).stats.tbl(2,2)];
    all_n_responsive = [all_n_responsive;parameter(p).stats.tbl(1,2)];
    
    all_perc_resistant = [all_perc_resistant;parameter(p).stats.tbl(2,2)/...
        (parameter(p).stats.tbl(2,2)+parameter(p).stats.tbl(2,1))];
    all_perc_responsive = [all_perc_responsive;parameter(p).stats.tbl(1,2)/...
        (parameter(p).stats.tbl(1,2)+parameter(p).stats.tbl(1,1))];
    
end

str_resistant = cell(length(all_n_resistant),1);
str_responsive = cell(length(all_n_resistant),1);

for i = 1:length(all_n_resistant)
    str_resistant{i} = sprintf('%d (%1.1f%%)',all_n_resistant(i),all_perc_resistant(i)*100);
    str_responsive{i} = sprintf('%d (%1.1f%%)',all_n_responsive(i),all_perc_responsive(i)*100);
end

eeg_supplemental_table = table(all_names,str_responsive,str_resistant,all_stat,all_p_value,...
    'VariableNames',{'Feature','Responsive','Resistant',...
    'OddsRatio','PValue'});

eeg_supplemental_table.Feature{1} = 'GSW awake';
eeg_supplemental_table.Feature{2} = 'GSW asleep';
eeg_supplemental_table.Feature{3} = 'PSW awake';
eeg_supplemental_table.Feature{4} = 'PSW asleep';
eeg_supplemental_table.Feature{5} = 'PST awake';
eeg_supplemental_table.Feature{6} = 'PST asleep';
eeg_supplemental_table.Feature{7} = 'GPFA awake';
eeg_supplemental_table.Feature{8} = 'GPFA asleep';
eeg_supplemental_table.Feature{9} = 'GLVFA awake';
eeg_supplemental_table.Feature{10} = 'GLVFA asleep';
eeg_supplemental_table.Feature{11} = 'Focal discharges awake';
eeg_supplemental_table.Feature{12} = 'Focal discharges asleep';
eeg_supplemental_table.Feature{13} = 'Focal slowing awake';
eeg_supplemental_table.Feature{14} = 'Focal slowing asleep';

eeg_supplemental_table.OddsRatio = arrayfun(@(x)sprintf('%1.2f',x),eeg_supplemental_table.OddsRatio,'UniformOutput',false);
eeg_supplemental_table.PValue = arrayfun(@(x)sprintf('%1.3f',x),eeg_supplemental_table.PValue,'UniformOutput',false);
writetable(eeg_supplemental_table,[results_folder,'SuppTable1.csv']);

%% Multivariate analysis
new_table.drug_resistant = logical(new_table.drug_resistant); % convert dr to logical
mdl = fitglm(new_table,'drug_resistant~gsw+psw+pst+gpfa+glvfa+foc_dis+foc_slow+duration_minutes',...
    'Distribution','binomial','CategoricalVars',{'gsw','psw','pst',...
    'gpfa','glvfa','foc_dis','foc_slow'});
mdl

% Simpler model
mdl2 = fitglm(new_table,'drug_resistant~pst+duration_minutes',...
    'Distribution','binomial','CategoricalVars',{'pst'});
mdl2


% Try taking the log the length of EEG
duration_log = log(new_table.duration_minutes);
new_table = addvars(new_table,duration_log);

mdl3 = fitglm(new_table,'drug_resistant~pst+duration_log',...
    'Distribution','binomial','CategoricalVars',{'pst'});
mdl3

%% Output variables to a csv file, to be read in R for FIRTH logistic regression using logistf
writetable(new_table,r_file_path);
    



end