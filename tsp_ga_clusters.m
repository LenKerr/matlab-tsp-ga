%TSP_GA_CLUSTERS Traveling Salesman Problem (TSP) Genetic Algorithm (GA) for Clusters
%   Finds a (near) optimal solution to the TSP by setting up a GA to search
%   for the shortest route (least distance for the salesman to travel to
%   each cluster exactly once and return to the starting cluster, while
%   visiting just a single unit per cluster)
%
% Summary:
%     1. A single salesman travels to each of the clusters and completes the
%        route by returning to the cluster he started from
%     2. Each cluster is visited by the salesman exactly once
%     3. Exactly one unit per cluster is visited by the salesman
%
% Input:
%     USERCONFIG (structure) with zero or more of the following fields:
%     - XY (float) is an Nx2 matrix of unit locations, where N is the number of cities
%     - DMAT (float) is an NxN matrix of point to point distances/finalCost
%     - CLUSTERS is an Nx1 cell array containing cluster->unit links (see example)
%     - POPSIZE (scalar integer) is the size of the population (should be divisible by 4)
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
%     - OPTROUTE (integer array) is the best route found by the algorithm (between clusters)
%     - OPTPATH (integer array) is the best route found by the algorithm (between units)
%     - MINDIST (scalar float) is the cost of the best route
%
% Usage:
%     tsp_ga_clusters
%       -or-
%     tsp_ga_clusters(userConfig)
%       -or-
%     resultStruct = tsp_ga_clusters;
%       -or-
%     resultStruct = tsp_ga_clusters(userConfig);
%       -or-
%     [...] = tsp_ga_clusters('Param1',Value1,'Param2',Value2, ...);
%
% Example:
%     % Let the function create an example problem to solve
%     tsp_ga_clusters;
%
% Example:
%     % Request the output structure from the solver
%     resultStruct = tsp_ga_clusters;
%
% Example:
%     nClusters = 15; nDims = 2;
%     unitsPerCluster = floor(5*rand(nClusters,1))+3;
%     nUnits = sum(unitsPerCluster);
%     xyCluster = 90*rand(nClusters,nDims);
%     clusters = cell(nClusters,1);
%     xy = zeros(nUnits,nDims);
%     c = 0; clrMat = hsv(nClusters);
%     figure; hold on;
%     for i = 1:nClusters
%         iUnits = unitsPerCluster(i);
%         index = c+(1:iUnits);
%         xyUnits = 10*rand(iUnits,nDims) + xyCluster(i(ones(iUnits,1)),:);
%         xy(index,:) = xyUnits;
%         clusters{i} = index;
%         c = c + iUnits;
%         plot(xyUnits(:,1),xyUnits(:,2),'.','Color',clrMat(i,:));
%     end
%     a = meshgrid(1:nUnits);
%     dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),nUnits,nUnits);
%     figure;imagesc(dmat);
%     userConfig = struct('xy',xy,'dmat',dmat,'clusters',{clusters});
%     resultStruct = tsp_ga_clusters(userConfig);
%
% See also: tsp_ga, tsp_nn, mtsp_ga
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
%
function varargout = tsp_ga_clusters(varargin)
    
    
    %
    % Initialize default configuration
    %
    defaultConfig.xy          = [];
    defaultConfig.dmat        = [];
    defaultConfig.clusters    = [];
    defaultConfig.popSize     = 20;
    defaultConfig.numIter     = 200;
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
    clusters    = configStruct.clusters;
    popSize     = configStruct.popSize;
    numIter     = configStruct.numIter;
    showProg    = configStruct.showProg;
    showStatus  = configStruct.showStatus;
    showResult  = configStruct.showResult;
    showWaitbar = configStruct.showWaitbar;
    if isempty(xy) || isempty(clusters)
        nClusters = 15; nDims = 2;
        unitsPerCluster = randi([3 7],nClusters,1);
        nUnits = sum(unitsPerCluster);
        xyCluster = 90*rand(nClusters,nDims);
        clusters = cell(nClusters,1);
        xy = zeros(nUnits,nDims);
        c = 0;
        for i = 1:nClusters
            iUnits = unitsPerCluster(i);
            index = c+(1:iUnits);
            xyUnits = 10*rand(iUnits,nDims) + xyCluster(i(ones(iUnits,1)),:);
            xy(index,:) = xyUnits;
            clusters{i} = index;
            c = c + iUnits;
        end
    end
    if isempty(dmat)
        nPoints = size(xy,1);
        a = meshgrid(1:nPoints);
        dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),nPoints,nPoints);
    end
    
    
    %
    % Process inputs
    %
    nClusters = length(clusters);
    unitsPerCluster = zeros(nClusters,1);
    for i = 1:nClusters
        unitsPerCluster(i) = length(clusters{i});
    end
    nUnits = sum(unitsPerCluster);
    [minUnits,minUnitIndex] = min(unitsPerCluster);
    unitStarts = clusters{minUnitIndex};
    
    
    %
    % Error checks
    %
    [N,dims] = size(xy);
    [nr,nc] = size(dmat);
    if (N ~= nr) || (N ~= nc) || (N ~= nUnits)
        error('??? Invalid XY, DMAT, or CLUSTER inputs')
    end
    n = nClusters;
    
    
    %
    % Sanity checks
    %
    popSize     = 4*ceil(popSize/4);
    numIter     = max(1,round(real(numIter(1))));
    showProg    = logical(showProg(1));
    showStatus  = logical(showStatus(1));
    showResult  = logical(showResult(1));
    showWaitbar = logical(showWaitbar(1));
    
    
    %
    % Initialize the population
    %
    pop = zeros(popSize,n);
    pop(1,:) = (1:n);
    for k = 2:popSize
        pop(k,:) = randperm(n);
    end
    
    
    %
    % Seed the algorithm with a previous result if available
    %
    if isfield(userConfig,'optRoute')
        optRoute = userConfig.optRoute;
        isValid = isequal(pop(1,:),sort(optRoute));
        if isValid
            pop(1,:) = optRoute;
        end
    end
    
    
    %
    % Select the colors for the clusters
    %
    pclr = ~get(0,'DefaultAxesColor');
    clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];
    if (nClusters > 5)
        clr = hsv(nClusters);
    end
    
    
    %
    % Run the GA
    %
    adjacencyHistory = cell(1,popSize);
    dijkstraHistory = cell(1,popSize);
    globalMin = Inf;
    totalDist = zeros(1,popSize);
    distHistory = NaN(1,numIter);
    tmpPop = zeros(4,n);
    newPop = zeros(popSize,n);
    [isClosed,isStopped,isCancelled] = deal(false);
    if showProg
        hFig = figure('Name','TSP_GA_CLUSTERS | Current Best Solution', ...
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
        % Evaluate each population member (calculate total distance)
        %
        for p = 1:popSize
            
            %
            % Build adjacency matrix from cluster info
            %
            A = false(N);
            A(clusters{pop(p,n)},clusters{pop(p,1)}) = true;
            for k = 2:n
                A(clusters{pop(p,k-1)},clusters{pop(p,k)}) = true;
            end
            dijkstraScore = Inf;
            for iUnit = 1:minUnits
                iStart = unitStarts(iUnit);
                [cost,path] = dijkstra_closed(A,dmat,iStart,iStart);
                if (cost < dijkstraScore)
                    dijkstraScore = cost;
                    dijkstraHistory{p} = path;
                end
            end
            adjacencyHistory{p} = A;
            totalDist(p) = dijkstraScore;
        end
        
        
        %
        % SELECT THE BEST
        %   This section of code finds the best solution in the current
        %   population and stores it if it is better than the previous best.
        %
        [minDist,index] = min(totalDist);
        distHistory(iter) = minDist;
        if (minDist < globalMin)
            globalMin = minDist;
            optRoute = pop(index,:);
            optPath = dijkstraHistory{index};
            if showProg
                
                %
                % Plot the best route
                %
                for i = 1:nClusters
                    c = clusters{i};
                    if (dims > 2), plot3(hAx,xy(c,1),xy(c,2),xy(c,3),'.','Color',clr(i,:));
                    else, plot(hAx,xy(c,1),xy(c,2),'.','Color',clr(i,:)); end
                    hold(hAx,'on');
                end
                rte = optPath;
                if (dims > 2), plot3(hAx,xy(rte,1),xy(rte,2),xy(rte,3),'o-','Color',pclr);
                else, plot(hAx,xy(rte,1),xy(rte,2),'o-','Color',pclr); end
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
        % Genetic algorithm operators
        %
        randomOrder = randperm(popSize);
        for p = 4:4:popSize
            rtes = pop(randomOrder(p-3:p),:);
            dists = totalDist(randomOrder(p-3:p));
            [ignore,idx] = min(dists); %#ok
            bestOf4Route = rtes(idx,:);
            routeInsertionPoints = sort(randperm(n,2));
            I = routeInsertionPoints(1);
            J = routeInsertionPoints(2);
            for k = 1:4 % Mutate the best to get three new routes
                tmpPop(k,:) = bestOf4Route;
                switch k
                    case 2 % Flip
                        tmpPop(k,I:J) = tmpPop(k,J:-1:I);
                    case 3 % Swap
                        tmpPop(k,[I J]) = tmpPop(k,[J I]);
                    case 4 % Slide
                        tmpPop(k,I:J) = tmpPop(k,[I+1:J I]);
                    otherwise % Do nothing
                end
            end
            newPop(p-3:p,:) = tmpPop;
        end
        pop = newPop;
        
        
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
    % Show the final results
    %
    if showResult
        
        %
        % Plot the GA results
        %
        figure('Name','TSP_GA_CLUSTERS | Results','Numbertitle','off');
        subplot(2,2,1);
        for i = 1:nClusters
            c = clusters{i};
            if (dims > 2), plot3(xy(c,1),xy(c,2),xy(c,3),'.','Color',clr(i,:));
            else, plot(xy(c,1),xy(c,2),'.','Color',clr(i,:)); end
            hold on
        end
        title('Unit Locations');
        subplot(2,2,2);
        clusterOrder = cat(2,clusters{optRoute});
        imagesc(dmat(clusterOrder,clusterOrder));
        title('Distance Matrix');
        subplot(2,2,3);
        for i = 1:nClusters
            c = clusters{i};
            if (dims > 2), plot3(xy(c,1),xy(c,2),xy(c,3),'.','Color',clr(i,:));
            else, plot(xy(c,1),xy(c,2),'.','Color',clr(i,:)); end
            hold on
        end
        rte = optPath;
        if (dims > 2), plot3(xy(rte,1),xy(rte,2),xy(rte,3),'o-','Color',pclr);
        else, plot(xy(rte,1),xy(rte,2),'o-','Color',pclr); end
        title(sprintf('Total Distance = %1.4f',minDist));
        subplot(2,2,4);
        plot(distHistory,'b','LineWidth',2);
        title('Best Solution History');
        set(gca,'XLim',[1 length(distHistory)],'YLim',[0 1.1*max([1 distHistory])]);
        
    end
    
    
    %
    % Return output
    %
    if nargout
        
        %
        % Create anonymous functions for plot generation
        %
        plotPoints  = @(s)plot(s.xy(:,1),s.xy(:,2),'k.');
        plotResult  = @(s)plot(s.xy(s.optSolution,1),s.xy(s.optSolution,2),'r.-');
        plotHistory = @(s)plot(s.distHistory,'b-','LineWidth',2);
        plotMatrix  = @(s)imagesc(s.dmat(s.optSolution,s.optSolution));
        
        
        %
        % Save results in output structure
        %
        resultStruct = struct( ...
            'xy',          xy, ...
            'dmat',        dmat, ...
            'clusters',    {clusters}, ...
            'popSize',     popSize, ...
            'numIter',     numIter, ...
            'showProg',    showProg, ...
            'showResult',  showResult, ...
            'showWaitbar', showWaitbar, ...
            'optRoute',    optRoute, ...
            'optPath',     optPath(1:nClusters), ...
            'optSolution', optPath([1:nClusters 1]), ...
            'plotPoints',  plotPoints, ...
            'plotResult',  plotResult, ...
            'plotHistory', plotHistory, ...
            'plotMatrix',  plotMatrix, ...
            'distHistory', distHistory, ...
            'minDist',     minDist);
        
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
            isRunning = false;
        else
            delete(hFig);
        end
    end
    
end

%
% Subfunction DIJKSTRA (Solves shortest path problem)
%
function [finalCost,finalPath] = dijkstra_closed(A,C,origin,target)
    
    
    %
    % Process inputs
    %
    n = size(A,1);
    [I,J] = find(A);
    E = [I J];
    
    
    %
    % Initializations
    %
    iTable = zeros(n,1);
    minCost = Inf(n,1);
    isSettled = false(n,1);
    paths = num2cell(NaN(n,1));
    I = origin;
    paths(I) = {I};
    
    
    %
    % Find the minimum cost and path using Dijkstra's Algorithm
    %
    while ~isSettled(target)
        
        %
        % Update the Table
        %
        jTable = iTable;
        iTable(I) = 0;
        neighborIds = find(E(:,1) == I);
        
        
        %
        % Calculate the costs to the neighbor points and record paths
        %
        for k = 1:length(neighborIds)
            J = E(neighborIds(k),2);
            if ~isSettled(J)
                cij = C(I,J);
                isNull = ~jTable(J);
                if isNull || (jTable(J) > (jTable(I) + cij))
                    iTable(J) = jTable(I) + cij;
                    paths{J} = [paths{I} J];
                else
                    iTable(J) = jTable(J);
                end
            end
        end
        
        K = find(iTable);
        
        %
        % Find the minimum value in the table
        %
        N = find(iTable(K) == min(iTable(K)));
        
        if isempty(N)
            break
        else
            
            %
            % Settle the minimum value
            %
            I = K(N(1));
            minCost(I) = iTable(I);
            isSettled(I) = true;
        end
        
    end
    
    
    %
    % Store costs and paths
    %
    finalCost = minCost(target);
    finalPath = paths{target};
    
end

