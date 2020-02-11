function parameter = ige_stats

%% Parameters

% simplify model
simple_model = 1; % if 1, this will remove focal features, slowing, and concatenate wake and sleep

% minimum number of elements in any cell of contingency table to do chi2
min_num_fisher = 10;

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

% tell it which parameters are the awake version of the one below it
eeg_parameters = [7:22,25:32];
awake_parameters = [7:2:21];
summary_parameters = [25:32];

% Select parameters to put in multivariate analysis
response = 1; % drug resistance is response variable (binary)
if simple_model == 1
    predictors = [3,24,25:29]; % sex, sw/psw/pst/gpfa/glvfa, duration
    pred_cat = [1,0,ones(1,length(predictors)-2)];
else
    predictors = [3,7:22,24]; % sex, all eeg features, duration
    pred_cat = [ones(1,length(predictors)-1),0]; % all categorical except duration
end


%% File path
csv_path = "/Users/erinconrad/Desktop/residency stuff/R25/ige project/data/IGEDatabase-DeidentifiedData_DATA_2020-02-10_1514.csv";

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
%{
sleep_col = find(strcmp(data.Properties.VariableNames,'sleep')); 
duration_col = find(strcmp(data.Properties.VariableNames,'eeg_duration')); 
exclude_eeg_col = find(strcmp(data.Properties.VariableNames,'exclude_eeg___0'));  
sex_col = find(strcmp(data.Properties.VariableNames,'sex'));
dr_col = find(strcmp(data.Properties.VariableNames,'drug_resistant'));
vpa_col = find(strcmp(data.Properties.VariableNames,'vpa'));
feature_cols = duration_col+1:duration_col+16; 
%}

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
    
    % now concatenate sleep and wake for new columns
    count = 0;
    for j = awake_parameters
        count = count + 1; % which of the awake parameters it is
        
        % get the columns corresponding to this AND the sleep version
        curr_cols = table2array(data(curr_eeg_rows,...
            parameter(j).column:parameter(j).column+1));
        
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


%% Identify exclusion rows
exclude = any([new_table.exclude_clinical___0==1,...
    new_table.exclude_eeg___0 == 1],2);

%% Univariate stats
% Get if drug resistant
dr_row = new_table.drug_resistant;
dr_row_logical = dr_row == 1;

% Loop through categorical parameters of interest
for p = [3:5,eeg_parameters]
    
    % get appropriate column
    name = parameter(p).name;
    
    % Get number who are resistant and number who have parameter
    resistant = new_table.drug_resistant;
    par = new_table.(name);
    
    % Find rows with nans or rows where resistant == 3
    nan_rows = any([isnan(resistant),isnan(par)],2);
    res_3_rows = resistant == 3;
    
    % summary of all unallowed rows
    all_unallowed = nan_rows|exclude|res_3_rows;
     
    
    % do chi squared
    [tbl,chi2,pval,labels] = crosstab(resistant(~all_unallowed),...
        par(~all_unallowed));
    
    % Don't do the analysis if all par is 0
    if length(unique(par(~all_unallowed))) == 1
        parameter(p).stats.tbl = tbl;
        parameter(p).stats.other = nan;
        continue
    end
    
    % if the value of any cell in the contingency table is <=min_num_fisher, do fisher
    % exact test instead
    if any(tbl(:)) <= min_num_fisher
        [~,pval,stats] = fishertest(tbl);
        parameter(p).stats.tbl = tbl;
        parameter(p).stats.p = pval;
        parameter(p).stats.OddsRatio = stats.OddsRatio;
        parameter(p).stats.ConfidenceInterval = stats.ConfidenceInterval;
    else
        parameter(p).stats.tbl = tbl;
        parameter(p).stats.chi2 = chi2;
        parameter(p).stats.p = pval;
        parameter(p).stats.labels = labels;
    end
    
    
    
    
end

% Wilcoxon rank sum for duration
p = 24;
dur = parameter(p).name;
[pval,~,stats] = ranksum(new_table.(dur)(dr_row_logical),new_table.(dur)(~dr_row_logical));
parameter(p).stats.p = pval;
parameter(p).stats.other = stats;
parameter(p).stats.avg_dur = [nanmean(new_table.(dur)(dr_row_logical)),nanmean(new_table.(dur)(~dr_row_logical))];

%% Multivariate analysis
% Get columns corresponding to response and predictor variables
response_col = parameter(response).column;
predictor_col = [];
for i = 1:length(predictors)
    predictor_col = [predictor_col,parameter(predictors(i)).column];
end

% Make new table with just the response and predictors, removing unallowed
% rows
tbl_for_mult = new_table(~all_unallowed,[response_col, predictor_col]);

% Convert things to arrays to make it friendly for the model
response_vec = table2array(tbl_for_mult(:,1));
predictor_array = table2array(tbl_for_mult(:,2:end));
cat_cols = pred_cat;

% Make response vector binary (right now it's 1s and 2s)
response_vec = logical(response_vec-1);


%% Output variables to a csv file, to be read in R for FIRTH logistic regression logistf

    
%% Multivariate analysis where we only include variables with p < 0.2 for univariate stats
% Find variables with p < 0.2
keep_predictors = zeros(length(predictors),1);
for i = 1:length(predictors)
    
    % if no p value, don't include it
    if isfield(parameter(predictors(i)).stats,'p') == 0
        continue;
    end
    
    % Get the pvalue
    pval = parameter(predictors(i)).stats.p;
    
    % Keep it if the p value is < 0.2
    if pval < 0.2
        keep_predictors(i) = 1;
    end
    
end

keep_predictors = logical(keep_predictors);


% Show the variable names
fprintf('\n\nPredictors included in the model:\n');
for i = 1:length(predictors)
    if keep_predictors(i) == 1
        fprintf('%s, p = %1.3f\n',parameter(predictors(i)).name,...
            parameter(predictors(i)).stats.p);
    end
end

% Do the model
mdl = fitglm(predictor_array(:,keep_predictors),response_vec,... % the first argument is the predictor, 2nd is the response
    'linear','Distribution','binomial','link','logit',... % establish that it's binomial
    'Categorical',logical(cat_cols(keep_predictors))); % tell it which predictors are categorical
mdl

mdl = fitglm(predictor_array(:,[2 5]),response_vec,... % the first argument is the predictor, 2nd is the response
    'linear','Distribution','binomial','link','logit',... % establish that it's binomial
    'Categorical',logical(cat_cols([2 5]))); 


end