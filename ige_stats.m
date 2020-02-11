function parameter = ige_stats

%% Parameters
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
parameter(24).name = 'duration_minutes';


%% File path
csv_path = "/Users/erinconrad/Desktop/residency stuff/R25/ige project/data/IGEDatabase-DeidentifiedData_DATA_2020-02-10_1514.csv";

%% Load csv file
data = readtable(csv_path);

%% Get column numbers
for i = 1:length(parameter)
    parameter(i).column = find(strcmp(data.Properties.VariableNames,parameter(i).name));
    if isempty(parameter(i).column) == 1
        if i == 24
            parameter(i).column = parameter(i-1).column;
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
    if curr_num == data.record_id(i)
        eeg_row(i) = 1;
    end
    curr_num = data.record_id(i);
end

% Prep array for duration in minutes
duration_minutes = zeros(size(eeg_row));

% Loop through rows and add eeg data to main clinical row
for i = 1:size(data,1)-1
    
    % continue if it's an eeg row
    if eeg_row(i) == 1, continue; end
    
    % Find the rows corresponding to the eeg data for this clinical row
    curr_eeg_rows = [];
    k = 1; % the possible eeg row we're checking
    while 1
        if i+k >size(data,1)
            break
        end
        if eeg_row(i+k) == 1
            curr_eeg_rows = [curr_eeg_rows;i+k];
        else
            break
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
        
        if length(dur_text) == 4
            dur_num = duration(dur_text,'InputFormat','hh:mm');
        elseif length(dur_text) == 8
            dur_num = duration(dur_text,'InputFormat','hh:mm:ss');
        end
        total_dur = dur_num + total_dur;

    end
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
for p = [3:5,7:22]
    
    % get appropriate column
    name = parameter(p).name;
    
    % Get number who are resistant and number who have parameter
    resistant = new_table.drug_resistant;
    par = new_table.(name);
    
    % Find rows with nans
    nan_rows = any([isnan(resistant),isnan(par)],2);
    
    % do chi squared
    [tbl,chi2,pval,labels] = crosstab(resistant(~nan_rows),par(~nan_rows));
    
    parameter(p).stats.tbl = tbl;
    parameter(p).stats.chi2 = chi2;
    parameter(p).stats.p = pval;
    parameter(p).stats.labels = labels;
    
    
end

% Wilcoxon rank sum for duration
p = 24;
dur = parameter(p).name;
[pval,~,stats] = ranksum(new_table.(dur)(dr_row_logical),new_table.(dur)(~dr_row_logical));
parameter(p).stats.pval = pval;
parameter(p).stats.other = stats;
parameter(p).stats.avg_dur = [nanmean(new_table.(dur)(dr_row_logical)),nanmean(new_table.(dur)(~dr_row_logical))];

end