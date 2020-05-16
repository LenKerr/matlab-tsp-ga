% PARAM_OVERRIDE
%
% Filename: get_config.m
%
% Description: Takes a default parameter structure along with user-defined
%     parameters and merges them by overriding the defaults with user values
%
% Date: 04/30/14
%
% Author:
%     Joseph Kirk
%     jdkirk630@gmail.com
%
% Inputs:
%     defaultConfig - structure containing default parameters
%     userConfig    - structure containing user-defined parameters
%
% Outputs:
%     config        - structure containing the full merged parameter set
%
% Usage:
%     config = get_config(defaultConfig,userConfig);
%
% See also: fieldnames
%
function config = get_config(defaultConfig,userConfig)
    
    
    %
    % Initialize the configuration structure as the default
    %
    config = defaultConfig;
    
    
    %
    % Extract the field names of the default configuration structure
    %
    defaultFields = fieldnames(defaultConfig);
    
    
    %
    % Extract the field names of the user configuration structure
    %
    userFields = fieldnames(userConfig);
    nUserFields = length(userFields);
    
    
    %
    % Override any default configuration fields with user values
    %
    for i = 1:nUserFields
        userField = userFields{i};
        isField = strcmpi(defaultFields,userField);
        if nnz(isField) == 1
            thisField = defaultFields{isField};
            config.(thisField) = userConfig.(userField);
        end
    end
    
end

