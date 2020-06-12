%TSP_NN Traveling Salesman Problem (TSP) Nearest Neighbor (NN) Algorithm
%   The Nearest Neighbor algorithm produces different results depending on
%   which city is selected as the starting point. This function determines
%   the Nearest Neighbor routes for multiple starting points and returns
%   the best of those routes
%
% Summary:
%     1. A single salesman travels to each of the cities and completes the
%        route by returning to the city he started from
%     2. Each city is visited by the salesman exactly once
%
% Input:
%     USERCONFIG (structure) with zero or more of the following fields:
%     - XY (float) is an Nx2 matrix of city locations, where N is the number of cities
%     - DMAT (float) is an NxN matrix of point to point distances/costs
%     - POPSIZE (scalar integer) is the size of the population (should be <= N)
%     - SHOWPROG (scalar logical) shows the GA progress if true
%     - SHOWRESULT (scalar logical) shows the GA results if true
%     - SHOWWAITBAR (scalar logical) shows a waitbar if true
%
% Input Notes:
%     1. Rather than passing in a structure containing these fields, any/all of
%        these inputs can be passed in as parameter/value pairs in any order instead.
%     2. Field/parameter names are case insensitive but must match exactly otherwise.
%
% Output:
%     RESULTSTRUCT (structure) with the following fields:
%         (in addition to a record of the algorithm configuration)
%     - OPTROUTE (integer array) is the best route found by the algorithm
%     - MINDIST (scalar float) is the cost of the best route
%
% Usage:
%     tsp_nn
%       -or-
%     tsp_nn(userConfig)
%       -or-
%     resultStruct = tsp_nn;
%       -or-
%     resultStruct = tsp_nn(userConfig);
%       -or-
%     [...] = tsp_nn('Param1',Value1,'Param2',Value2, ...);
%
% Example:
%     % Let the function create an example problem to solve
%     tsp_nn;
%
% Example:
%     % Request the output structure from the solver
%     resultStruct = tsp_nn;
%
% Example:
%     % Pass a random set of user-defined XY points to the solver
%     userConfig = struct('xy',10*rand(50,2));
%     resultStruct = tsp_nn(userConfig);
%
% Example:
%     % Pass a more interesting set of XY points to the solver
%     n = 100;
%     phi = (sqrt(5)-1)/2;
%     theta = 2*pi*phi*(0:n-1);
%     rho = (1:n).^phi;
%     [x,y] = pol2cart(theta(:),rho(:));
%     xy = 10*([x y]-min([x;y]))/(max([x;y])-min([x;y]));
%     userConfig = struct('xy',xy);
%     resultStruct = tsp_nn(userConfig);
%
% Example:
%     % Pass a random set of 3D (XYZ) points to the solver
%     xyz = 10*rand(50,3);
%     userConfig = struct('xy',xyz);
%     resultStruct = tsp_nn(userConfig);
%
% Example:
%     % Turn off the plots but show a waitbar
%     userConfig = struct('showProg',false,'showResult',false,'showWaitbar',true);
%     resultStruct = tsp_nn(userConfig);
%
% See also: tsp_ga, tspo_ga, tspof_ga, tspofs_ga
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
%
function varargout = tsp_nn(varargin)
    
    
    %
    % Initialize default configuration
    %
    defaultConfig.xy          = 10*rand(100,2);
    defaultConfig.dmat        = [];
    defaultConfig.popSize     = Inf;
    defaultConfig.showProg    = true;
    defaultConfig.showStatus  = true;
    defaultConfig.showResult  = true;
    defaultConfig.showWaitbar = false;
    
    
    %
    % Interpret user configuration inputs
    %
    if ~nargin
        userConfig = struct();
    elseif isstruct(varargin{1})
        userConfig = varargin{1};
    else
        try
            userConfig = struct(varargin{:});
        catch
            error('??? Expected inputs are either a structure or parameter/value pairs');
        end
    end
    
    
    %
    % Override default configuration with user inputs
    %
    configStruct = get_config(defaultConfig,userConfig);
    
    
    %
    % Extract configuration
    %
    xy          = configStruct.xy;
    dmat        = configStruct.dmat;
    popSize     = configStruct.popSize;
    showProg    = configStruct.showProg;
    showStatus  = configStruct.showStatus;
    showResult  = configStruct.showResult;
    showWaitbar = configStruct.showWaitbar;
    if isempty(dmat)
        nPoints = size(xy,1);
        a = meshgrid(1:nPoints);
        dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),nPoints,nPoints);
    end
    
    
    %
    % Verify inputs
    %
    [N,dims] = size(xy);
    [nr,nc] = size(dmat);
    if (N ~= nr) || (N ~= nc)
        error('??? Invalid XY or DMAT inputs')
    end
    n = N;
    
    
    %
    % Sanity checks
    %
    popSize     = max(1,min(n,round(real(popSize(1)))));
    showProg    = logical(showProg(1));
    showStatus  = logical(showStatus(1));
    showResult  = logical(showResult(1));
    showWaitbar = logical(showWaitbar(1));
    
    
    %
    % Initialize the population
    %
    pop = zeros(popSize,n);
    
    
    %
    % Run the NN
    %
    distHistory = NaN(1,popSize);
    [isClosed,isStopped,isCancelled] = deal(false);
    if showProg
        hFig = figure('Name','TSP_NN | Current Solution', ...
            'Numbertitle','off','CloseRequestFcn',@close_request);
        hAx = gca;
        if showStatus
            [hStatus,isCancelled] = figstatus(0,popSize,[],hFig);
        end
    end
    if showWaitbar
        hWait = waitbar(0,'Searching for near-optimal solution ...', ...
            'CreateCancelBtn',@cancel_search);
    end
    isRunning = true;
    for p = 1:popSize
        d = 0;
        thisRte = zeros(1,n);
        isVisited = false(1,n);
        I = p;
        isVisited(I) = true;
        thisRte(1) = I;
        for k = 2:n
            dists = dmat(I,:);
            dists(logical(isVisited)) = NaN;
            dMin = min(dists(~isVisited));
            J = find(dists == dMin,1);
            isVisited(J) = true;
            thisRte(k) = J;
            d = d + dmat(I,J);
            I = J;
        end
        d = d + dmat(I,p);
        pop(p,:) = thisRte;
        distHistory(p) = d;
        
            
        %
        % Plot the current route
        %
        if showProg
            rte = thisRte([1:n 1]);
            if (dims > 2), plot3(hAx,xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
            else, plot(hAx,xy(rte,1),xy(rte,2),'r.-'); end
            title(hAx,sprintf('Total Distance = %1.4f, Iteration = %d', ...
                distHistory(p),p));
            drawnow;
        end
        
        
        %
        % Update the status bar and check cancellation status
        %
        if showProg && showStatus && ~mod(p,ceil(popSize/100))
            [hStatus,isCancelled] = figstatus(p,popSize,hStatus,hFig);
        end
        if (isStopped || isCancelled)
            break
        end
        
        
        %
        % Update the waitbar
        %
        if showWaitbar && ~mod(p,ceil(popSize/325))
            waitbar(p/popSize,hWait);
        end
        
    end
    if showProg && showStatus
        figstatus(popSize,popSize,hStatus,hFig);
    end
    if showWaitbar
        delete(hWait);
    end
    isRunning = false;
    if isClosed
        delete(hFig);
    end
    
    
    %
    % Find the minimum distance route
    %
    [minDist,index] = min(distHistory);
    optRoute = pop(index,:);
    
    
    %
    % Format the optimal solution
    %
    index = find(optRoute == 1,1);
    optSolution = [optRoute([index:n 1:index-1]) 1];
    
    
    %
    % Show the final results
    %
    if showResult
        
        %
        % Plot the best route
        %
        if showProg && ~isClosed
            rte = optRoute([1:n 1]);
            if (dims > 2), plot3(hAx,xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
            else, plot(hAx,xy(rte,1),xy(rte,2),'r.-'); end
            title(hAx,sprintf('Total Distance = %1.4f',minDist));
        end
        
        %
        % Plot the NN results
        %
        figure('Name','TSP_NN | Results','Numbertitle','off');
        subplot(2,2,1);
        pclr = ~get(0,'DefaultAxesColor');
        if (dims > 2), plot3(xy(:,1),xy(:,2),xy(:,3),'.','Color',pclr);
        else, plot(xy(:,1),xy(:,2),'.','Color',pclr); end
        title('City Locations');
        subplot(2,2,2);
        imagesc(dmat(optRoute,optRoute));
        title('Distance Matrix');
        subplot(2,2,3);
        rte = optRoute([1:n 1]);
        if (dims > 2), plot3(xy(rte,1),xy(rte,2),xy(rte,3),'r.-');
        else, plot(xy(rte,1),xy(rte,2),'r.-'); end
        title(sprintf('Total Distance = %1.4f',minDist));
        subplot(2,2,4);
        plot(sort(distHistory,'descend'),'b','LineWidth',2);
        title('Distances');
        set(gca,'XLim',[0 popSize+1],'YLim',[0 1.1*max([1 distHistory])]);
    end
    
    
    %
    % Return output
    %
    if nargout
        
        %
        % Create anonymous functions for plot generation
        %
        plotPoints  = @(s)plot(s.xy(:,1),s.xy(:,2),'.','Color',~get(gca,'Color'));
        plotResult  = @(s)plot(s.xy(s.optSolution,1),s.xy(s.optSolution,2),'r.-');
        plotHistory = @(s)plot(s.distHistory,'b-','LineWidth',2);
        plotMatrix  = @(s)imagesc(s.dmat(s.optSolution,s.optSolution));
        
        
        %
        % Save results in output structure
        %
        resultStruct = struct( ...
            'xy',          xy, ...
            'dmat',        dmat, ...
            'popSize',     popSize, ...
            'showProg',    showProg, ...
            'showResult',  showResult, ...
            'showWaitbar', showWaitbar, ...
            'optRoute',    optRoute, ...
            'optSolution', optSolution, ...
            'plotPoints',  plotPoints, ...
            'plotResult',  plotResult, ...
            'plotHistory', plotHistory, ...
            'plotMatrix',  plotMatrix, ...
            'distHistory', distHistory, ...
            'minDist',     minDist, ...
            'pop',         pop);
        
        varargout = {resultStruct};
        
    end
    
    
    %
    % Nested function to cancel search
    %
    function cancel_search(varargin)
        isStopped = true;
    end
    
    
    %
    % Nested function to close the figure window
    %
    function close_request(varargin)
        if isRunning
            [isClosed,isStopped] = deal(true);
        else
            delete(hFig);
        end
    end
    
end

