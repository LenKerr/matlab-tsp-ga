% MTSP_RS Multiple Traveling Salesmen Problem (M-TSP) Random Search (RS) Algorithm
%   The Random Search algorithm generates random solutions and evaluates
%   their fitness
%
% Summary:
%     1. Each salesman travels to a unique set of cities and completes the
%        route by returning to the city he started from
%     2. Each city is visited by exactly one salesman
%
% Input:
%     USERCONFIG (structure) with zero or more of the following fields:
%     - XY (float) is an Nx2 matrix of city locations, where N is the number of cities
%     - DMAT (float) is an NxN matrix of city-to-city distances or costs
%     - NSALESMEN (scalar integer) is the number of salesmen to visit the cities
%     - MINTOUR (scalar integer) is the minimum tour length for any of the salesmen
%     - POPSIZE (scalar integer) is the size of the population (should be <= N)
%     - NUMITER (scalar integer) is the number of desired iterations for the algorithm to run
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
%     - OPTBREAK (integer array) is the list of route break points (these specify the indices
%         into the route used to obtain the individual salesman routes)
%     - MINDIST (scalar float) is the total distance traveled by the salesmen
%
% Route/Breakpoint Details:
%     If there are 10 cities and 3 salesmen, a possible route/break
%     combination might be: rte = [5 6 9 1 4 2 8 10 3 7], brks = [3 7]
%     Taken together, these represent the solution [5 6 9][1 4 2 8][10 3 7],
%     which designates the routes for the 3 salesmen as follows:
%         . Salesman 1 travels from city 5 to 6 to 9 and back to 5
%         . Salesman 2 travels from city 1 to 4 to 2 to 8 and back to 1
%         . Salesman 3 travels from city 10 to 3 to 7 and back to 10
%
% Usage:
%     mtsp_rs
%       -or-
%     mtsp_rs(userConfig)
%       -or-
%     resultStruct = mtsp_rs;
%       -or-
%     resultStruct = mtsp_rs(userConfig);
%       -or-
%     [...] = mtsp_rs('Param1',Value1,'Param2',Value2, ...);
%
% Example:
%     % Let the function create an example problem to solve
%     mtsp_rs;
%
% Example:
%     % Request the output structure from the solver
%     resultStruct = mtsp_rs;
%
% Example:
%     % Pass a random set of user-defined XY points to the solver
%     userConfig = struct('xy',10*rand(35,2));
%     resultStruct = mtsp_rs(userConfig);
%
% Example:
%     % Pass a more interesting set of XY points to the solver
%     n = 50;
%     phi = (sqrt(5)-1)/2;
%     theta = 2*pi*phi*(0:n-1);
%     rho = (1:n).^phi;
%     [x,y] = pol2cart(theta(:),rho(:));
%     xy = 10*([x y]-min([x;y]))/(max([x;y])-min([x;y]));
%     userConfig = struct('xy',xy);
%     resultStruct = mtsp_rs(userConfig);
%
% Example:
%     % Pass a random set of 3D (XYZ) points to the solver
%     xyz = 10*rand(35,3);
%     userConfig = struct('xy',xyz);
%     resultStruct = mtsp_rs(userConfig);
%
% Example:
%     % Turn off the plots but show a waitbar
%     userConfig = struct('showProg',false,'showResult',false,'showWaitbar',true);
%     resultStruct = mtsp_rs(userConfig);
%
% See also: mtsp_ga, tsp_ga, tsp_rs
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
%
function varargout = mtsp_rs(varargin)
    
    
    %
    % Initialize default configuration
    %
    defaultConfig.xy          = 10*rand(30,2);
    defaultConfig.dmat        = [];
    defaultConfig.nSalesmen   = 5;
    defaultConfig.minTour     = 3;
    defaultConfig.popSize     = Inf;
    defaultConfig.numIter     = 5e3;
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
    nSalesmen   = configStruct.nSalesmen;
    minTour     = configStruct.minTour;
    popSize     = configStruct.popSize;
    numIter     = configStruct.numIter;
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
    nSalesmen   = max(1,min(n,round(real(nSalesmen(1)))));
    minTour     = max(1,min(floor(n/nSalesmen),round(real(minTour(1)))));
    popSize     = max(1,min(n,round(real(popSize(1)))));
    numIter     = max(1,round(real(numIter(1))));
    showProg    = logical(showProg(1));
    showStatus  = logical(showStatus(1));
    showResult  = logical(showResult(1));
    showWaitbar = logical(showWaitbar(1));
    
    
    %
    % Initializations for route break point selection
    %
    nBreaks = nSalesmen-1;
    dof = n - minTour*nSalesmen;          % degrees of freedom
    addto = ones(1,dof+1);
    for k = 2:nBreaks
        addto = cumsum(addto);
    end
    cumProb = cumsum(addto)/sum(addto);
    
    
    %
    % Initialize the populations
    %
    popRoute = zeros(popSize,n);         % population of routes
    popBreak = zeros(popSize,nBreaks);   % population of breaks
    popRoute(1,:) = (1:n);
    popBreak(1,:) = rand_breaks();
    for k = 2:popSize
        popRoute(k,:) = randperm(n);
        popBreak(k,:) = rand_breaks();
    end
    
    
    %
    % Seed the algorithm with a previous result if available
    %
    if all(isfield(userConfig,{'optRoute','optBreak'}))
        optRoute = userConfig.optRoute;
        optBreak = userConfig.optBreak;
        isValidRoute = isequal(popRoute(1,:),sort(optRoute));
        isValidBreak = all(optBreak > 0) && all(optBreak <= n) && ...
            (length(optBreak) == nBreaks) && ~any(mod(optBreak,1));
        if isValidRoute && isValidBreak
            popRoute(1,:) = optRoute;
            popBreak(1,:) = optBreak;
        end
    end
    
    
    %
    % Select the colors for the plotted routes
    %
    pclr = ~get(0,'DefaultAxesColor');
    clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];
    if (nSalesmen > 5)
        clr = hsv(nSalesmen);
    end
    
    
    %
    % Run the Random Search (RS)
    %
    globalMin = Inf;
    distHistory = NaN(1,numIter);
    [isClosed,isStopped,isCancelled] = deal(false);
    if showProg
        hFig = figure('Name','MTSP_RS | Current Best Solution', ...
            'Numbertitle','off','CloseRequestFcn',@close_request);
        hAx = gca;
        if showStatus
            [hStatus,isCancelled] = figstatus(0,numIter,[],hFig);
        end
    end
    if showWaitbar
        hWait = waitbar(0,'Searching for near-optimal solution ...', ...
            'CreateCancelBtn',@cancel_search);
    end
    isRunning = true;
    for iter = 1:numIter
        
        %
        % EVALUATE SOLUTIONS
        %   This section of code computes the total cost of each solution
        %   in the population. The actual code that gets executed uses a
        %   much faster (vectorized) method to calculate the route lengths
        %   compared to the triple for-loop below (provided for reference)
        %   but gives the same result.
        %
        %     totalDist = zeros(popSize,1);
        %     for p = 1:popSize
        %         d = 0;
        %         pRoute = popRoute(p,:);
        %         pBreak = popBreak(p,:);
        %         rng = [[1 pBreak+1];[pBreak n]]';
        %         for s = 1:nSalesmen
        %             d = d + dmat(pRoute(rng(s,2)),pRoute(rng(s,1)));
        %             for k = rng(s,1):rng(s,2)-1
        %                 d = d + dmat(pRoute(k),pRoute(k+1));
        %             end
        %         end
        %         totalDist(p) = d;
        %     end
        %
        row = popRoute;
        col = popRoute(:,[2:n 1]);
        for p = 1:popSize
            brk = popBreak(p,:);
            col(p,[brk n]) = row(p,[1 brk+1]);
        end
        ind = N*(col-1) + row;
        totalDist = sum(dmat(ind),2);
        
        
        %
        % SELECT THE BEST
        %   This section of code finds the best solution in the current
        %   population and stores it if it is better than the previous best.
        %
        [minDist,index] = min(totalDist);
        distHistory(iter) = minDist;
        if (minDist < globalMin)
            globalMin = minDist;
            optRoute = popRoute(index,:);
            optBreak = popBreak(index,:);
            rng = [[1 optBreak+1];[optBreak n]]';
            if showProg
                
                %
                % Plot the best route
                %
                for s = 1:nSalesmen
                    rte = optRoute([rng(s,1):rng(s,2) rng(s,1)]);
                    if (dims > 2), plot3(hAx,xy(rte,1),xy(rte,2),xy(rte,3),'.-','Color',clr(s,:));
                    else, plot(hAx,xy(rte,1),xy(rte,2),'.-','Color',clr(s,:)); end
                    hold(hAx,'on');
                end
                title(hAx,sprintf('Total Distance = %1.4f, Iteration = %d',minDist,iter));
                hold(hAx,'off');
                drawnow;
            end
        end
        
        
        %
        % Update the status bar and check cancellation status
        %
        if showProg && showStatus && ~mod(iter,ceil(numIter/100))
            [hStatus,isCancelled] = figstatus(iter,numIter,hStatus,hFig);
        end
        if (isStopped || isCancelled)
            break
        end
        
        
        %
        % MODIFY THE POPULATION
        %   This section of code randomly searches for new solutions
        %
        popRoute(1,:) = optRoute;
        popBreak(1,:) = optBreak;
        for k = 2:popSize
            popRoute(k,:) = randperm(n);
            popBreak(k,:) = rand_breaks();
        end
        
        
        %
        % Update the waitbar
        %
        if showWaitbar && ~mod(iter,ceil(numIter/325))
            waitbar(iter/numIter,hWait);
        end
        
    end
    if showProg && showStatus
        figstatus(numIter,numIter,hStatus,hFig);
    end
    if showWaitbar
        delete(hWait);
    end
    isRunning = false;
    if isClosed
        delete(hFig);
    end
    
    
    %
    % Append prior distance history if present
    %
    if isfield(userConfig,'distHistory')
        priorHistory = userConfig.distHistory;
        isNan = isnan(priorHistory);
        distHistory = [priorHistory(~isNan) distHistory];
    end
    
    
    %
    % Format the optimal solution
    %
    optSolution = cell(nSalesmen,1);
    rng = [[1 optBreak+1];[optBreak n]]';
    for s = 1:nSalesmen
        optSolution{s} = optRoute([rng(s,1):rng(s,2) rng(s,1)]);
    end
    
    
    %
    % Show the final results
    %
    if showResult
        
        %
        % Plot the GA results
        %
        figure('Name','MTSP_RS | Results','Numbertitle','off');
        subplot(2,2,1);
        if (dims > 2), plot3(xy(:,1),xy(:,2),xy(:,3),'.','Color',pclr);
        else, plot(xy(:,1),xy(:,2),'.','Color',pclr); end
        title('City Locations');
        subplot(2,2,2);
        imagesc(dmat(optRoute,optRoute));
        title('Distance Matrix');
        subplot(2,2,3);
        for s = 1:nSalesmen
            rte = optSolution{s};
            if (dims > 2), plot3(xy(rte,1),xy(rte,2),xy(rte,3),'.-','Color',clr(s,:));
            else, plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:)); end
            title(sprintf('Total Distance = %1.4f',minDist));
            hold on;
        end
        subplot(2,2,4);
        plot(distHistory,'b','LineWidth',2);
        title('Best Solution History');
        set(gca,'YLim',[0 1.1*max([1 distHistory])]);
    end
    
    
    %
    % Return output
    %
    if nargout
        
        %
        % Create anonymous functions for plot generation
        %
        plotPoints  = @(s)plot(s.xy(:,1),s.xy(:,2),'.','Color',~get(gca,'Color'));
        plotResult  = @(s)cellfun(@(s,i)plot(s.xy(i,1),s.xy(i,2),'.-', ...
            'Color',rand(1,3)),repmat({s},size(s.optSolution)),s.optSolution);
        plotHistory = @(s)plot(s.distHistory,'b-','LineWidth',2);
        plotMatrix  = @(s)imagesc(s.dmat(cat(2,s.optSolution{:}),cat(2,s.optSolution{:})));
        
        
        %
        % Save results in output structure
        %
        resultStruct = struct( ...
            'xy',          xy, ...
            'dmat',        dmat, ...
            'nSalesmen',   nSalesmen, ...
            'minTour',     minTour, ...
            'popSize',     popSize, ...
            'numIter',     numIter, ...
            'showProg',    showProg, ...
            'showResult',  showResult, ...
            'showWaitbar', showWaitbar, ...
            'optRoute',    optRoute, ...
            'optBreak',    optBreak, ...
            'optSolution', {optSolution}, ...
            'plotPoints',  plotPoints, ...
            'plotResult',  plotResult, ...
            'plotHistory', plotHistory, ...
            'plotMatrix',  plotMatrix, ...
            'distHistory', distHistory, ...
            'minDist',     minDist);
        
        varargout = {resultStruct};
        
    end
    
    
    %
    % Generate random set of break points
    %
    function breaks = rand_breaks()
        if (minTour == 1) % No Constraints on breaks
            tmpBreaks = randperm(n-1);
            breaks = sort(tmpBreaks(1:nBreaks));
        else % Force breaks to be at least the minimum tour length
            nAdjust = find(rand < cumProb,1)-1;
            spaces = randi(nBreaks,1,nAdjust);
            adjust = zeros(1,nBreaks);
            for kk = 1:nBreaks
                adjust(kk) = sum(spaces == kk);
            end
            breaks = minTour*(1:nBreaks) + cumsum(adjust);
        end
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

