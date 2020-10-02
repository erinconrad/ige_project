function mdl = feature_selection(tbl,resp,vars)

%{
This function selects the features for a multivariate logistic regression
model to predict drug resistance. I starts with those variables that had
p<0.1 on univariate analysis for association with drug resistance. I first
do a univariate logistic regression with those variables and sort the
variables by model AIC. I then stepwise add variables in ascending order of
AIC, and if the multivariate model has a lower AIC I add the variable and
if it doesn't I reject the variable.
%}

%% Get AIC for each variable taken alone
aic_uni = zeros(length(vars),1);
for i = 1:length(aic_uni)
    mdl = fitglm(tbl,'ResponseVar',resp,...
        'PredictorVars',vars{i});
    aic = mdl.ModelCriterion.AIC;
    aic_uni(i) = aic;
end

%% Sort them in order of their AIC (lowest to highest)
[sorted_aics,I] = sort(aic_uni);
vars = vars(I);

%% Go through and add variables one by one and reject if AIC goes up
% In order of lowest to highest AIC
curr_aic = sorted_aics(1); % start with lowest AIC
curr_vars = vars(1);
for i = 2:length(aic_uni)
    curr_vars{end+1} = vars{i}; % Add variable to current list of variables
    mdl = fitglm(tbl,'ResponseVar',resp,...
        'PredictorVars',curr_vars);
    aic = mdl.ModelCriterion.AIC;
    if aic < curr_aic
        curr_aic = aic; % change lowest aic to this if it is indeed lower
    else
        curr_vars(end) = []; % if not, reject this variable
    end    
end

%% Do model on these final features
mdl = fitglm(tbl,'ResponseVar',resp,...
        'PredictorVars',curr_vars);
    

end