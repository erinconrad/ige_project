%{
This function performs the analysis for the paper "Using generalized 
polyspike train to predict drug-resistant idiopathic generalized epilepsy".

It, and de-identified data needed to run it, can be found at:
https://github.com/erinconrad/ige_project

The survival analysis is done in R on code that can be found in the same
github repository.

Erin Conrad, 2020
%}

function ige_stats

%% Parameters for users to change
doing_from_github = 0; % Set to 1 if running the code from the github repository (most people)

% File paths
csv_path = 'IGEDatabase-DeidentifiedData_DATA_2020-04-10_1654.csv'; % point to the data
results_folder = '../results/'; % point to where you want to save results
r_data_path = '../data/'; % where to save a table of data to be used for R analysis

display_random_entries = 1; % Set to 1 to display data for random patients, to double check that it agrees with data table

%% Other parameters

% Add path to extra tools
if doing_from_github == 0
    addpath(genpath('../tools/'));
end

% Get U and H statistic? (These employ versions of built-in Matlab code for
% the Wilcoxon rank sum test that I modified to give
% the U and H statistics. I could not commit this to github because the
% software is proprietary from Matlab, and so if running this from the
% github repository this should be zero).
get_u_and_h = ~doing_from_github;

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
parameter(34).name = 'gsw_time';
parameter(35).name = 'psw_time';
parameter(36).name = 'pst_time';
parameter(37).name = 'gpfa_time';
parameter(38).name = 'total_time_first_gsw'; 
parameter(39).name = 'total_time_first_psw';
parameter(40).name = 'total_time_first_pst';
parameter(41).name = 'total_time_first_gpfa';
parameter(42).name = 'total_time_first_pst_or_gpfa';
parameter(43).name = 'pst_or_gpfa';
parameter(44).name = 'pst_or_gpfa_awake';
parameter(45).name = 'pst_or_gpfa_asleep';


% Pretty names for figure
parameter(38).pretty_name = 'GSW'; 
parameter(39).pretty_name = 'PSW';
parameter(40).pretty_name = 'PST';
parameter(41).pretty_name = 'GPFA';

% tell it which parameters are the awake version of the one below it
eeg_parameters = [7:22,25:32];
awake_parameters = 7:2:21;
summary_parameters = 25:32;

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
fprintf('\nPrepping data table...');
new_table = data;

total_eegs = 0;

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

% Prep arrays for time to first occurrence of various eeg features
first_feature = nan(length(eeg_row),4);

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
            total_eegs = total_eegs + 1;
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
            %fprintf('Warning, duration text for row %d says %s\n',i,dur_text);
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
        
        % Find time to first occurrence of various eeg features
        if ismember(count,1:4) % loop through features we're checking (in order, gsw, psw, pst, gpfa)
            
            % Skip it if they don't have the feature
            if any(curr_cols(:)) == 0
                continue
            end
            
            time_elapsed = 0;
            found_feature_time = 0;
            
            % Get order of EEG rows by date of EEG (age at eeg so de-id)
            eeg_dates = table2array(data(curr_eeg_rows,...
            parameter(33).column)); %33 is age_at_eeg
            [~,date_order] = sort(eeg_dates);
            
            % Loop through EEG rows by date of EEG
            for d = curr_eeg_rows(date_order)'
                
                % Get the time to feature for this eeg
                eeg_time = table2array(data(d,parameter(34-1+count).column)); %34 is gsw time
                eeg_time = duration(eeg_time{1},'InputFormat','hh:mm');
                
                % If we find the feature in this EEG, fill array and break
                if ~isnan(eeg_time)
                    first_feature(i,count) = (time_elapsed) + minutes(eeg_time);
                    found_feature_time = 1;
                    break
                else
                % If we don't find it, add duration of the EEG
                    cur_dur = table2array(data(d,parameter(2).column));
                    cur_dur = cur_dur{1};
        
                    % Fix for mistakes in text
                    cur_dur = strrep(cur_dur,'.',':');

                    if length(cur_dur) == 4 || length(cur_dur) == 5
                        cur_dur = duration(cur_dur,'InputFormat','hh:mm');
                    elseif length(cur_dur) == 8
                        cur_dur = duration(cur_dur,'InputFormat','hh:mm:ss');
                    elseif isempty(cur_dur) == 1
                        cur_dur = 0;
                    elseif strcmp(cur_dur,'41 minutes') == 1 % correction for mistake in text
                        cur_dur = duration('00:41','InputFormat','hh:mm');
                    else
                        %fprintf('Warning, duration text for row %d says %s\n',i,cur_dur);
                        cur_dur = 0;
                    end
                    if cur_dur == 0
                        time_elapsed = time_elapsed + (cur_dur);
                    else
                        time_elapsed = time_elapsed + minutes(cur_dur);
                    end
                
                end
            end
            % If we haven't found the time after this loop and not excluding, throw an error
            if found_feature_time == 0
                if any(table2array(data(curr_eeg_rows,parameter(6).column))) == 0 && ...
                        table2array(data(min(curr_eeg_rows)-1,parameter(23).column)) == 0
                    error('what\n'); 
                end
            end
            
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

% Add time to first occurrence info
total_time_first_gsw = first_feature(:,1);
total_time_first_psw = first_feature(:,2);
total_time_first_pst = first_feature(:,3);
total_time_first_gpfa = first_feature(:,4);
total_time_first_pst_or_gpfa = min(first_feature(:,3:4),[],2);
new_table = addvars(new_table,total_time_first_gsw);
new_table = addvars(new_table,total_time_first_psw);
new_table = addvars(new_table,total_time_first_pst);
new_table = addvars(new_table,total_time_first_gpfa);
new_table = addvars(new_table,total_time_first_pst_or_gpfa);

% Add gpt or gpfa
pst_or_gpfa = any([new_table.pst,new_table.gpfa],2);
pst_or_gpfa_awake = any([new_table.pst___1,new_table.gpfa___1],2);
pst_or_gpfa_asleep = any([new_table.pst___2,new_table.gpfa___2],2);

new_table = addvars(new_table,pst_or_gpfa);
new_table = addvars(new_table,pst_or_gpfa_awake);
new_table = addvars(new_table,pst_or_gpfa_asleep);

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
    error('\nWarning, record ID %d with zero duration. Removing...\n',...
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
fprintf('Done.\n----------------------------------------------------\n');

fprintf('\n\nNow displaying random entries to test that it agrees with original data table:\n\n');
if display_random_entries
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
fprintf('\n----------------------------------------------------\n');
fprintf('Now doing statistics...\n\n');
%% Relationship between PST and VPA (I don't report this)
% do chi squared
[tbl,chi2,pval,labels] = crosstab(new_table.pst,...
    new_table.vpa);

% Do Fisher exact test
[~,pval,stats] = fishertest(tbl);


%% Wilcoxon rank sums for duration and age at first eeg
fprintf('Comparing continuous parameters between drug resistant and drug responsive patients:\n\n');
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
for p = [25:32,43]
    % get appropriate column
    name = parameter(p).name;
    
    % Get number who have parameter
    par = new_table.(name);
    
    if sum(par) == 0, continue; end
    
    [pval,~,stats] = ranksum(new_table.duration_minutes(par == 1),...
    new_table.duration_minutes(par == 0));
    if get_u_and_h == 1
        [alt_pval,~,U] = ranksum_output_u(new_table.duration_minutes(par == 1),...
        new_table.duration_minutes(par == 0));
    else
        alt_pval = pval;
        U = nan;
    end

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
duration_feature_table.Feature{8} = 'PST or GPFA';
duration_feature_table.AverageDurationWithFeature = ...
    arrayfun(@(x)sprintf('%1.1f',x),duration_feature_table.AverageDurationWithFeature,'UniformOutput',false);
duration_feature_table.AverageDurationWithoutFeature = ...
    arrayfun(@(x)sprintf('%1.1f',x),duration_feature_table.AverageDurationWithoutFeature,'UniformOutput',false);

% Add asterisks to p
p_stars = cell(length(duration_feature_table.p),1);
for i = 1:length(duration_feature_table.p)
    if duration_feature_table.p(i) < 0.001/length(duration_feature_table.p)
        p_stars{i} = '***';
    elseif duration_feature_table.p(i) < 0.01/length(duration_feature_table.p)
        p_stars{i} = '**';
    elseif duration_feature_table.p(i) < 0.05/length(duration_feature_table.p)
        p_stars{i} = '*'; 
    else
        p_stars{i} = ''; 
    end
end

duration_feature_table.p = arrayfun(@(x)sprintf('%1.3f',x),duration_feature_table.p,'UniformOutput',false);
for i = 1:length(duration_feature_table.p)
    if str2num(duration_feature_table.p{i}) < 0.001
        duration_feature_table.p{i} = ['<0.001',p_stars{i}];
    else
        duration_feature_table.p{i} = [duration_feature_table.p{i},p_stars{i}];
    end
end


%% Statistical test comparing time to first occurrence of feature
fprintf('\n----------------------------------------------------\n');
fprintf('\n\nComparing time to first occurrence of different features:\n\n');
% Get times to first occurrences
all_times = [new_table.total_time_first_gsw,new_table.total_time_first_psw,...
    new_table.total_time_first_pst,new_table.total_time_first_gpfa];


% Kruskall wallis with post-hoc comparisons
[p,tbl,stats] = kruskalwallis(all_times,[],'off');

% Post-hoc dunn-sidak tests
c = multcompare(stats,'CType','dunn-sidak','Display','off');

% Print the stats
fprintf(['The median time to first occurrence of feature is:\n'...
    '%1.1f minutes for GSW,\n'...
    '%1.1f minutes for PSW,\n'...
    '%1.1f minutes for PST,\n'...
    '%1.1f minutes for GPFA.\n'...
    'Kruskall-Wallis p = %1.3f\n'],nanmedian(new_table.total_time_first_gsw),...
    nanmedian(new_table.total_time_first_psw),...
    nanmedian(new_table.total_time_first_pst),...
    nanmedian(new_table.total_time_first_gpfa),...
    p);    

%% Figure for percentage of patients with feature by certain times
% X axis is duration of eeg
% Y axis is number of patients whose eeg is that duration or lower who have
% the feature of interest
fprintf('\n----------------------------------------------------\n');
fprintf('\n\nMaking figure showing percentage of patients with feature by certain times...\n\n');
fig1 = figure;
set(fig1,'position',[100 378 1324 422])

subplot(1,2,1)
legend_names = {};

features_to_plot = 38:41; % this is time to first feature for gsw, psw, pst, gpfa

dur_short = zeros(length(features_to_plot),1);
count = 0;
for p = features_to_plot
    count = count + 1;
    
    % get appropriate column
    name = parameter(p).name;
    pretty_name = parameter(p).pretty_name;
    legend_names = [legend_names,pretty_name];
    
    % Get first time at which that feature occurs
    par = new_table.(name);
    
    % Define duration step size
    dur_steps = 0:20:max(par)+20;
    n_par_dur = zeros(1,length(dur_steps));
    perc_par_dur = zeros(1,length(dur_steps));
    
    % Loop through duration steps
    for d = 1:length(dur_steps)
        
        % Get the number of patients who have parameter by
        % less than or equal to specified duration
        n_par_dur(d) = sum(par<=dur_steps(d));
        
        % Get the percentage of patients who have parameter within that
        % specified duration (divide by the total number with that
        % parameter)
        perc_par_dur(d) = sum(par<=dur_steps(d))/sum(~isnan(par))*100;
    end
    
    plot(dur_steps/60,perc_par_dur,'linewidth',2)
    hold on
    
    % Compare less than 1 hour EEGs
    dur_short(count) = sum(par<=1*60)/sum(~isnan(par))*100;
    
end

xlim([0 24])
xlabel('Time to first feature occurrence (hours)')
ylabel({'Percentage'});
legend(legend_names,'location','southeast')
set(gca,'fontsize',20)
annotation('textbox',[0.07 0.88 0.1 0.1],'String','A',...
    'linestyle','none','fontsize',35);


subplot(1,2,2)
bar(dur_short)
hold on
xticklabels(legend_names)
set(gca,'fontsize',20)
xlabel('EEG feature')
ylabel({'Percentage';'with occurrence in <1 hour'});
annotation('textbox',[0.48 0.88 0.1 0.1],'String','B',...
    'linestyle','none','fontsize',35);
print(fig1,[results_folder,'Figure1'],'-depsc')



%% Categorical Univariate stats
% Loop through categorical parameters of interest
for p = [3:5,eeg_parameters,43:45]
    
    % Note that male == 1, female == 0
    % Tried VPA == 1, not tried VPA == 0
    
    % get appropriate column
    name = parameter(p).name;
    
    % Get number who have parameter
    par = new_table.(name);
    
    % do chi squared to get table
    [tbl,chi2,pval,labels] = crosstab(new_table.drug_resistant,...
        par);
    
    % Don't do the analysis if all par is 0
    if length(unique(par)) == 1
        parameter(p).stats.tbl = [tbl,[0;0]];
        parameter(p).stats.other = nan;
        continue
    end
    
    % Fisher exact test (the one I actually use given small numbers)
    [~,pval,stats] = fishertest(tbl);
    parameter(p).stats.labels = labels;
    parameter(p).stats.tbl = tbl;
    parameter(p).stats.p = pval;
    parameter(p).stats.OddsRatio = stats.OddsRatio;
    parameter(p).stats.ConfidenceInterval = stats.ConfidenceInterval;
    
 
end

%% Clinical table

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

fprintf('\n----------------------------------------------------\n');
fprintf('\n\nClinical table:\n');
clinical_table
fprintf('\n');
writetable(clinical_table,[results_folder,'Table1.csv']);

fprintf('\n----------------------------------------------------\n');
fprintf('\n\nFeature duration table:\n');
duration_feature_table
fprintf('\n');
writetable(duration_feature_table,[results_folder,'Table2.csv']);

%% Make a summary table of univariate statistics for eeg features
all_names = {};
all_p_value = [];
all_stat = [];
all_n_resistant = [];
all_n_responsive = [];
all_ci = [];

all_perc_resistant = [];
all_perc_responsive = [];

for p = [25:32,43]
    if strcmp(parameter(p).name,'duration_minutes') == 1, continue; end
    if strcmp(parameter(p).name,'sex') == 1, continue; end
    
    if isfield(parameter(p).stats,'p') == 0
        continue
    end
    
    all_names = [all_names;parameter(p).name];
    all_p_value = [all_p_value;parameter(p).stats.p];
    all_stat = [all_stat;parameter(p).stats.OddsRatio];
    all_ci = [all_ci;parameter(p).stats.ConfidenceInterval];

    
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

% Convert the odds ratio to a string with CI
or_str = cell(length(all_stat),1);
for i = 1:length(or_str)
    if isnan(all_stat(i))
        or_str{i} = sprintf('NaN');
    elseif all_stat(i) == inf
        or_str{i} = sprintf('Inf');
    else
        or_str{i} = sprintf('%1.1f (%1.1f-%1.1f)',all_stat(i),all_ci(i,1),all_ci(i,2));
    end
end

eeg_feature_table = table(all_names,str_responsive,str_resistant,or_str,all_p_value,...
    'VariableNames',{'Feature','Responsive','Resistant',...
    'OddsRatio','PValue'});

% Format table
%eeg_feature_table.OddsRatio = arrayfun(@(x)sprintf('%1.2f',x),eeg_feature_table.OddsRatio,'UniformOutput',false);
eeg_feature_table.PValue = arrayfun(@(x)sprintf('%1.3f',x),eeg_feature_table.PValue,'UniformOutput',false);
eeg_feature_table.Feature{1} = 'GSW';
eeg_feature_table.Feature{2} = 'PSW';
eeg_feature_table.Feature{3} = 'PST';
eeg_feature_table.Feature{4} = 'GPFA';
eeg_feature_table.Feature{5} = 'GLVFA';
eeg_feature_table.Feature{6} = 'Focal discharges';
eeg_feature_table.Feature{7} = 'Focal slowing';
eeg_feature_table.Feature{8} = 'PST or GPFA';
fprintf('\n----------------------------------------------------\n');
fprintf('\n\nEEG feature table:\n');
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
all_ci = [];

for p = [7:16,19:22,44:45] %wake and sleep eeg features
    if strcmp(parameter(p).name,'duration_minutes') == 1, continue; end
    if strcmp(parameter(p).name,'sex') == 1, continue; end
    
    
    
    all_names = [all_names;parameter(p).name];
    if isfield(parameter(p).stats,'p') == 0
        all_p_value = [all_p_value;nan];
        all_stat = [all_stat;nan];
        all_ci = [all_ci;nan nan];
    else
        all_p_value = [all_p_value;parameter(p).stats.p];
        all_stat = [all_stat;parameter(p).stats.OddsRatio];
        all_ci = [all_ci;parameter(p).stats.ConfidenceInterval];
    end
    

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

% Convert the odds ratio to a string with CI
or_str = cell(length(all_stat),1);
for i = 1:length(or_str)
    if isnan(all_stat(i))
        or_str{i} = sprintf('NaN');
    elseif all_stat(i) == inf
        or_str{i} = sprintf('Inf');
    else
        or_str{i} = sprintf('%1.1f (%1.1f-%1.1f)',all_stat(i),all_ci(i,1),all_ci(i,2));
    end
end

eeg_supplemental_table = table(all_names,str_responsive,str_resistant,or_str,all_p_value,...
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
eeg_supplemental_table.Feature{15} = 'PST or GPFA awake';
eeg_supplemental_table.Feature{16} = 'PST or GPFA asleep';


eeg_supplemental_table.PValue = arrayfun(@(x)sprintf('%1.3f',x),eeg_supplemental_table.PValue,'UniformOutput',false);
writetable(eeg_supplemental_table,[results_folder,'SuppTable1.csv']);


%% Is PST more frequent in sleep than wakefulness?
if 0
fprintf('%d patients (%1.1f%%) had PST while awake, and %d (%1.1f%%) had PST while asleep.\n',...
    sum(new_table.pst___1),sum(new_table.pst___1)/length(new_table.pst___1)*100,...
    sum(new_table.pst___2),sum(new_table.pst___2)/length(new_table.pst___2)*100);
end

%% Stratify by duration (<1 hour, 1-24 hours, >=24 hours)
fprintf('\n----------------------------------------------------\n');
fprintf('\n\nNow controlling for duration by stratification:\n');
sub_hour = new_table.duration_minutes <= 60;

[tbl] = crosstab(new_table.drug_resistant(sub_hour),...
    new_table.pst(sub_hour));
[~,pval,stats] = fishertest(tbl);

fprintf(['\nLooking only at sub-hour EEGs,\n'...
    '%d drug-responsive patients (%1.1f%%) have PST and %d resistant patients (%1.1f%%) have PST\n'...
    '(Fisher exact test: odds-ratio %1.1f, p = %1.3f)\n'],tbl(1,2),tbl(1,2)/(tbl(1,2)+tbl(1,1))*100,...
    tbl(2,2),tbl(2,2)/(tbl(2,2)+tbl(2,1))*100,stats.OddsRatio,pval);

[tbl] = crosstab(new_table.drug_resistant(~sub_hour),...
    new_table.pst(~sub_hour));
[~,pval,stats] = fishertest(tbl);

fprintf(['\nLooking only at super-hour EEGs,\n'...
    '%d drug-responsive patients (%1.1f%%) have PST and %d resistant patients (%1.1f%%) have PST\n'...
    '(Fisher exact test: odds-ratio %1.1f, p = %1.3f)\n'],tbl(1,2),tbl(1,2)/(tbl(1,2)+tbl(1,1))*100,...
    tbl(2,2),tbl(2,2)/(tbl(2,2)+tbl(2,1))*100,stats.OddsRatio,pval);


%% Log-rank test for survival analysis
fprintf('\n----------------------------------------------------\n');
fprintf('\n\nNow doing survival analysis:\n');
% survival analysis: analyzes the expected duration of time until one or more events happen
% The log-rank test compares the survival times of two or more groups
% Censoring is a form of missing data problem in which time to event is not observed

% I will do a log-rank test to compare the time to first PST, and consider
% the end of the EEG as censoring
x = nan(size(new_table,1),2);
x(:,1) = new_table.total_time_first_pst; % time to first pst (nan if never happens)
x(:,2) = isnan(new_table.total_time_first_pst); % 1 if no PST before end of EEG
x(isnan(new_table.total_time_first_pst),1) = new_table.duration_minutes(isnan(new_table.total_time_first_pst)); % set the time for the nan rows to be the EEG duration

% separate drug resistant and drug responsive
x1 = x(new_table.drug_resistant==1,:);
x2 = x(new_table.drug_resistant==0,:);

% My Kaplan-Meier plot
fprintf('Making plot...\n')
log_rank_erin(x1,x2,results_folder,'PST');

% prep for R (log rank test done in R)
all_x = nan(size(new_table,1),3);

% concatenate drug resistant and drug responsive together
all_x(:,1) = [x1(:,1);x2(:,1)]; 
all_x(:,2) = [x1(:,2);x2(:,2)]; 
all_x(:,3) = [ones(size(x1,1),1);zeros(size(x2,1),1)]; % the third column indicates if drug resistant

% This table has the survival time (either EEG duration or time to first
% PST) in the first column, whether PST was observed in the second column,
% and drug resistance in the third column
r_table = table(all_x(:,1),~all_x(:,2),all_x(:,3),'VariableNames',{'survtime','observed','resistant'});

% export table for R
writetable(r_table,[r_data_path,'r_table_pst.csv']);

fprintf('\n----------------------------------------------------\n');
fprintf('End of file. All tables and figures should be saved to your results folder.\nThe log-rank test can be done in R using the data file saved to your data folder.\n\n');

end