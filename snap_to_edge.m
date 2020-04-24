function [newpoints,Gmag,Gdir] = snap_to_edge(img,points,varargin)
%SNAP_TO_EDGE Snaps points to nearest edge of 2D image I.
% img     : 2D image
% points  : n x 2 set of points or a cell array with n x 2 point sets.
% Optional inputs (as 'argument',<value> pairs)
% spacing : pixel spacing in mm. Default: [1 1]
% maxd    : maximum distance to move points (in mm). Default: 1
% sigma   : standard deviation of the gaussian image filter. Default: 1
% method  : interpolation method. Default: 'linear'.
% direction: direction in which points can be moved: 'normal' (in direction
%            normal to curve) or 'gradient' (orthogonal to the image gradient)
%
% Bart Bolsterlee
% Neuroscience Research Australia
% 17/04/2020

%% Parse inputs
p = inputParser;
addRequired(p,'img')
addRequired(p,'points')
addParameter(p,'maxd',1,@isscalar)
addParameter(p,'spacing',[1 1])
addParameter(p,'sigma',1,@(x) isscalar(x) || isempty(x))
addParameter(p,'method','linear')
addParameter(p,'direction','gradient')
parse(p,img,points,varargin{:});

method  = p.Results.method; % interpolation method
spacing = p.Results.spacing;
maxd    = p.Results.maxd;
sigma   = p.Results.sigma;
direction = p.Results.direction;

PlotFlag = false; % Set to true for diagnostics purposes only.
%% Crop image
% Build up x and y-vector of image.
imdim = size(img);
xvec = (0 : imdim(2)-1) * spacing(2);
yvec = (0 : imdim(1)-1) * spacing(1);

flag = false;
if isnumeric(points)
    % Convert points to cell.
    flag = true;
    points = {points};
end

% Crop image to location where points are. This will speed up image
% gradient calculation.
% Note: this will give a problem if points are near the edge of the image...
padding = 3*maxd;
minx = min(cellfun(@(x) min(x(:,1)),points));
maxx = max(cellfun(@(x) max(x(:,1)),points));
miny = min(cellfun(@(x) min(x(:,2)),points));
maxy = max(cellfun(@(x) max(x(:,2)),points));
sel_x = xvec > (minx - padding) & xvec < (maxx + padding);
sel_y = yvec > (miny - padding) & yvec < (maxy + padding);


%% Calculate gradient direction and magnitude
if isempty(sigma)
    Ifilt = img(sel_y,sel_x);
else
    Ifilt = imgaussfilt(img(sel_y,sel_x),sigma);
end    
[Gmag,Gdir] = imgradient(Ifilt);

%%

n = 1001; % sample n points along the direction vector
d = linspace(-maxd,maxd,n)';
newpoints = cell(size(points));
if PlotFlag == true
    figure('Color','w')
    subplot(1,2,1)
    xlabel('Distance from original location (mm)')
    ylabel('Gradient magnitude')
end
for i = 1 : numel(points)
    % Calculate edge vectors in the gradient direction.
    N = NaN(size(points{i}));
    switch direction
        case 'gradient'
            N(:,1) = interp2(xvec(sel_x),yvec(sel_y),-cosd(Gdir),points{i}(:,1),points{i}(:,2),method);
            N(:,2) = interp2(xvec(sel_x),yvec(sel_y), sind(Gdir),points{i}(:,1),points{i}(:,2),method);
        case 'normal'
            [~,~,~,~,N] = fit_spline(points{i});
    end
    
    Nx = d*N(:,1)' + ones(n,1)*points{i}(:,1)';
    Ny = d*N(:,2)' + ones(n,1)*points{i}(:,2)';
    
    % Interpolate image at points along vector
    G_n = interp2(xvec(sel_x),yvec(sel_y),Gmag,Nx,Ny,method);
    
    % Find all peaks of curve and use the one closest to 0 (i.e. closest to
    % the position before snapping. This prevents snapping to the other
    % side of thin objects.
    nPoints = size(G_n,2);
    max_value = NaN(nPoints,1);
    idx       = NaN(nPoints,1);
    for j = 1 : nPoints
        % Index of all maxima of the curve
        max_idx = find(islocalmax(G_n(:,j))); 
        if isempty(max_idx)
            % No edge found within maximum distance. Keep the original
            % position.
            idx(j) = find(d==0);
        else
            % Get values of all maxima of the curve.
            all_max_value = G_n(max_idx,j);
            
            % Exclude local maxima that are less than half as high as the
            % global maximum. This should prevent snapping to trivial
            % edges.
            excl = all_max_value < 0.5*max(all_max_value);
            max_idx(excl) = [];
            [~,closest_idx] = min(abs(d(max_idx))); % index of closes maxima
            idx(j) = max_idx(closest_idx);
        end
        max_value(j) = G_n(idx(j),j);

    end
    if PlotFlag == true
        hold on
        plot(d,G_n,'k-')
        plot(d(idx),max_value,...
            'o',...
            'MarkerSize',8,...
            'MarkerFaceColor','r',...
            'MarkerEdgeColor','k')
    end
    
    % Calculate new position
    newpoints{i} = points{i} + d(idx)*[1 1] .* N;    
end

% Plot for diagnostics:
if PlotFlag == true
    subplot(1,2,2)
    image(xvec(sel_x),yvec(sel_y),Gmag)
    colormap gray
    hold on
    for i = 1 : length(points)
        h1 = plot(points{i}(:,1),points{i}(:,2),'-oy',...
            'LineWidth',2,...
            'MarkerFaceColor','y',...
            'MarkerEdgeColor','none',...
            'MarkerSize',6);
        h2 = plot(newpoints{i}(:,1),newpoints{i}(:,2),'-or',...
            'LineWidth',2,...
            'MarkerFaceColor','none',...
            'MarkerEdgeColor','r',...
            'MarkerSize',7);
    end
    legend([h1 h2],{'before','after'});
    axis equal off tight
end
if flag == true
    newpoints = newpoints{1};
end


end

