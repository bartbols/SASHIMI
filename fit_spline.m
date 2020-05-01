function [Ps,Ns,ts,t,N] = fit_spline( pos,varargin )
%FIT_SPLINE uses csape to fit a closed curve around the points given
%by P (n x 2 array).
%
% Optional input arguments as 'argument',<value> pairs:
% n:      Number of points on the spline curve . At these
%          points, the normal vectors are calculated as well.
%       
% conds:  'periodic' or 'clamped'. End conditions for slopes. Use
%         'periodic' to fit a closed curve. Default: periodic
% 
%
% OUTPUT:
% Vs : location of points sampled along the fitted curve
% Ns : normal vectors of points sampled along the fitted curve
% t  : t-values of control points
% ts : t-values of points sampled along the curve
% N  : normal vectors at control point location
%
% Bart Bolsterlee
% Neuroscience Research Australia
% April 2020

%% Parse inputs
p = inputParser;
addRequired(p,'pos')
addParameter(p,'n',[])
addParameter(p,'conds','periodic',@(x) any(strcmp(x,{'clamped','periodic'})))
parse(p,pos,varargin{:});
conds = p.Results.conds;
n     = p.Results.n;
%%

% The first and last point should be the same. If this is not the case, add
% the first point to the end of the array.
flag = false;
if ~isequal(pos(1,:),pos(end,:)) && strcmp(conds,'periodic')
    flag = true;
    pos = [pos;pos(1,:)];
end

% t-values of control points
t = cumsum(sqrt([0,diff(pos(:,1)')].^2 + [0,diff(pos(:,2)')].^2));
if isempty(n)
    % define number of samples based on the length of the curve so that
    % each segment is approximately 0.5 mm long.
    n = 2 * ceil(max(t));
end
% conds = 'variational';
px = csape(t,pos(:,1),conds);
py = csape(t,pos(:,2),conds);

% t-values of sampled points along the curve
ts = linspace(0,max(t),n);
Ps(:,1) = fnval(px,ts);
Ps(:,2) = fnval(py,ts);

if nargout > 1
    % Normal vector at sampled points.
    Ns(:,1) = -fnval(fnder(py),linspace(0,max(t),n));
    Ns(:,2) =  fnval(fnder(px),linspace(0,max(t),n));
    
    % Normalize to unit vector.
    Ns = Ns ./ (sqrt(sum(Ns.^2,2))*[1 1]);
    if nargout > 4
        % Also calculate the normal vectors at the original point
        % locations.
        N(:,1) = -fnval(fnder(py),t);
        N(:,2) =  fnval(fnder(px),t);
        
        % normalize to unit vector
        N = N ./ (sqrt(sum(N.^2,2))*[1 1]);
        if flag == true
            % Remove the last point again so that it matches with the
            % input.
            N = N(1:end-1,:);
            t = t(1:end-1);
        end
    end
end

end

