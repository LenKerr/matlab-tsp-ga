%FIGSTATUS Creates a status bar in a figure to show loop progress
%
% Description: Displays loop progress as a status bar in a figure window
%
% Inputs:
%     i       - (numeric) current loop index
%     n       - (numeric) total number of loop iterations
%     hStatus - (optional) handle to status object
%     hFig    - (optional) handle to figure
%
% Outputs:
%     hStatus     - handle to status object
%     isCancelled - scalar logical indicating user has clicked "cancel" button
%
% Usage:
%     [hStatus,isCancelled] = figstatus(i,n);
%       -or-
%     hStatus = figstatus(i,n);
%       -or-
%     figstatus(i,n,hStatus);
%       -or-
%     figstatus(i,n,hStatus,hFig);
%
% Example:
%     n = 1e3;
%     hStatus = figstatus(0,n);
%     for i = 1:n
%         % Computations go here ...
%         figstatus(i,n,hStatus);
%     end
%
% Example:
%     n = 1e7;
%     hStatus = figstatus(0,n);
%     for i = 1:n
%         % Computations go here ...
%         if ~mod(i,ceil(n/100))
%             figstatus(i,n,hStatus);
%         end
%     end
%
% Example:
%     n = 1e7;
%     figstatus(0,n);
%     for i = 1:n
%         % Computations go here ...
%         if ~mod(i,ceil(n/100))
%             figstatus(i,n);
%         end
%     end
%
% Example:
%     % Create a button to quit early if desired
%     n = 1e3;
%     [hStatus,isCancelled] = figstatus(0,n);
%     for i = 1:n
%         % Computations go here ...
%         [hStatus,isCancelled] = figstatus(i,n);
%         if isCancelled
%             % Normal use in a function would be to break out of the loop, but
%             %   if this example is run in the Command Window, just throw error
%             % break
%             error('??? User killed processing');
%         end
%     end
%
% See also: progress, uipanel, uicontrol, patch
%
function varargout = figstatus(i,n,hStatus,hFig)
    
    
    %
    % Check for between 2 and 4 inputs
    %
    narginchk(2,4);
    
    
    %
    % Set persistent quit trigger
    %
    persistent IS_CANCELLED;
    
    
    %
    % Initialize cancel button trigger
    %
    if isempty(IS_CANCELLED)
        IS_CANCELLED = false;
    end
    
    
    %
    % Send quit signal if triggered
    %
    if IS_CANCELLED
        IS_CANCELLED = false;
        varargout = {[],true};
        return
    end
    
    
    %
    % Check for status handle
    %
    if (nargin < 3) || isempty(hStatus) || ~ishandle(hStatus)
        
        %
        % Get current figure if not provided
        %
        if (nargin < 4) || isempty(hFig) || ~ishandle(hFig)
            hFig = gcf;
        end
        
        %
        % Look for handle to status
        %
        hStatus = findobj(hFig,'Tag','_figstatus_');
        if isempty(hStatus) && (i ~= n)
            
            %
            % Create status panel and bar
            %
            figColor = get(hFig,'Color');
            hStatusPanel = uipanel(hFig, ...
                'Units',           'normalized', ...
                'Position',        [0 0 1 0.05], ...
                'BackgroundColor', figColor, ...
                'Tag',             '_figstatuspanel_', ...
                'BorderType',      'etchedin');
            hStatusAx = axes( ...
                'Color',           'w', ...
                'Units',           'normalized', ...
                'Position',        [0 0 1 1], ...
                'XLim',            [0 1], ...
                'YLim',            [0 1], ...
                'Visible',         'off', ...
                'Parent',          hStatusPanel);
            if ~verLessThan('matlab','9.5') % R2018b
                set(hStatusAx.Toolbar,'Visible','off');
                % disableDefaultInteractivity(hStatusAx);
            end
            if (nargout > 1)
                set(hStatusAx,'Position',[0 0 0.95 1]);
                uicontrol(hStatusPanel, ...
                    'Units',       'normalized', ...
                    'Position',    [0.95 0 0.05 1], ...
                    'String',      'X', ...
                    'Callback',    @trigger_cancel);
            end
            % cmw = 0.1 * repmat([7 10 10 7],3,1)';
            % patch([0 0 1 1],[0 1 1 0],[0.5 0.5 0.5], ...
            %     'FaceVertexCData',cmw, ...
            %     'FaceColor','interp', ...
            %     'EdgeColor','none', ...
            %     'Parent',hStatusAx);
            grn = 0.1 * [2 6 2; 8 10 8];
            cmg = grn([1 2 2 1],:);
            hStatus = patch([0 0 1 1],[0 1 1 0],[0.5 0.5 0.5], ...
                'FaceVertexCData',cmg, ...
                'FaceColor','interp', ...
                'EdgeColor','none', ...
                'Tag','_figstatus_', ...
                'Parent',hStatusAx);
            
        end
    end
    
    
    %
    % Set status bar length
    %
    frac = (i / n);
    set(hStatus,'XData',[0 0 frac frac]);
    drawnow();
    
    
    %
    % If complete, cleanup by deleting entire status panel
    %
    if (i == n)
        hStatusPanel = findobj(hFig,'Tag','_figstatuspanel_');
        delete(hStatusPanel);
    end
    
    
    %
    % Pass output handle if requested
    %
    if nargout
        varargout = {hStatus,IS_CANCELLED};
    end
    
    
    %
    % Nested function to trigger quit
    %
    function trigger_cancel(varargin)
        IS_CANCELLED = true;
        hStatusPanel  = findobj(hFig,'Tag','_figstatuspanel_');
        delete(hStatusPanel);
    end
    
end

