function h = fig(varargin)
    if mod(numel(varargin),2)
        f = varargin{1};
        parameters = varargin(2:end);
    else
        parameters = varargin;
        f = [];
    end
    
    if ishandle(f)
        set(0,'CurrentFigure',f);
    elseif ~isempty(f)
        f = figure(f);
    else
        f = figure; shg;
    end
    
    if ~isempty(parameters)
        set(f,parameters{:});
    end
    
    if nargout == 1
        h = f;
    end
end