function varargout = register_slices(img1,img2,varargin)
%REGISTER_SLICES Registers two 2D images and then applies the estimated
%transformation field to transform the label of image 1 to image 2. This
%function uses Matlab's built-in demons registration (see doc imregdemons
%for details). The first input is the moving image and the second input
%the fixed image. The label should be defined on the moving image.
%
% Bart Bolsterlee, Neuroscience Research Australia
% March 2020

%% Parse inputs
p = inputParser;
addRequired(p,'img1')
addRequired(p,'img2')
% addRequired(p,'label1')
addParameter(p,'label1',[])
addParameter(p,'points',[])
% Registration settings
addParameter(p,'N',100,@isscalar)
addParameter(p,'AccumulatedFieldSmoothing',1.5)
addParameter(p,'PyramidLevels',3)
addParameter(p,'DisplayWaitbar',false)
addParameter(p,'spacing',[1 1])
addParameter(p,'method','linear')
addParameter(p,'crop',false)
addParameter(p,'padding',5,@isscalar)
parse(p,img1,img2,varargin{:});

% label1 = p.Results.label1;
N       = p.Results.N; % number of iterations
% demons settings. See 'doc imregdemons' for details.
afs     = p.Results.AccumulatedFieldSmoothing;
pl      = p.Results.PyramidLevels;
dw      = p.Results.DisplayWaitbar;
crop    = p.Results.crop;
padding = p.Results.padding;
label1  = p.Results.label1; % label defined on the moving image
points  = p.Results.points; % points in fixed image coordinates
spacing = p.Results.spacing; % pixel spacing in x and y direction
method  = p.Results.method; % interpolation method

label1_tf = [];
points_tf = [];
%% Register images

% first input  : moving image
% second input : fixed image (reference)
% D            : displacement vectors from the fixed image grid to a
%                corresponding location in the moving image in image pixel
%                dimensions.
% img1_reg     : warped moving image (should match the fixed image)

imdim = size(img2);
xvec = (0 : imdim(2)-1) * spacing(2);
yvec = (0 : imdim(1)-1) * spacing(1);

if crop == true && ~isempty(points)
    % Register only regions in which points are
    if iscell(points)
        minx = min(cellfun(@(x) min(x(:,1)),points));
        maxx = max(cellfun(@(x) max(x(:,1)),points));
        miny = min(cellfun(@(x) min(x(:,2)),points));
        maxy = max(cellfun(@(x) max(x(:,2)),points));
    else
        minx = min(points(:,1));
        maxx = max(points(:,1));
        miny = min(points(:,2));
        maxy = max(points(:,2));
    end
    sel_x = xvec > (minx - padding) & xvec < (maxx + padding);
    sel_y = yvec > (miny - padding) & yvec < (maxy + padding);
else
    sel_x = true(size(xvec));
    sel_y = true(size(yvec));

end
img1_reg = NaN(size(img2));
% D        = NaN([size(img2),2]);
[D,I_reg] = imregdemons(img1(sel_y,sel_x),img2(sel_y,sel_x),N,...
    'AccumulatedFieldSmoothing',afs,...
    'PyramidLevels',pl,...
    'DisplayWaitbar',dw);
img1_reg(sel_y,sel_x,:) = I_reg;
% D(sel_y,sel_x,:) = D2;

if ~isempty(label1)
    % Estimate the transformed label. label1_tf is label 1 transformed to
    % match image 2.
    label1_tf = imwarp(label1,D,'nearest');
end

% Transform the points
if ~isempty(points)
    if iscell(points)
        % Cell array with multiple point sets is provided. Transform all 
        % points per cell.
        points_tf = cell(size(points));
        for i = 1 : numel(points)
            if ~isempty(points{i})
                points_tf{i} = get_new_location(points{i});
            end
        end
    else
        % One array with points is provided.        
        points_tf = get_new_location(points);
    end
end

% Sub-functions
    function pnew = get_new_location(pts)        
            d = zeros(size(pts));
            d(:,1) = interp2(xvec(sel_x),yvec(sel_y),D(:,:,1),pts(:,1),pts(:,2),method); % x displacement in pixel coordinates
            d(:,2) = interp2(xvec(sel_x),yvec(sel_y),D(:,:,2),pts(:,1),pts(:,2),method); % y displacement in pixel coordinates
            % Update the new position vectors. The displacement field
            % preditect by imregdemons is in pixel coordinates, so it needs
            % to be converted back to image coordinates (flipped x/y axis,
            % and multiplied by pixel size).
            pnew = pts + d.*spacing([2 1]);        
    end
%%
if nargout > 0
    varargout{1} = img1_reg;
    if nargout > 1
        varargout{2} = label1_tf;
        if nargout > 2
            varargout{3} = points_tf;
        end
    end
end


end

