function varargout = register_slices(img1,img2,varargin)
%REGISTER_SLICES Registers two 2D images and then applies the estimated
%transformation field to transform the points of image 1 to image 2. This
%function uses Matlab's built-in demons registration (see doc imregdemons
%for details). The first input is the moving image and the second input
%the fixed image.
%
% INPUT
% ----- Required inputs -----
% img1    : reference ('fixed') image
% img2    : moving image
%
% ----- Optional inputs -----
% Optional input arguments are provided as 'argument',<value> pairs:
% points: nx2 array (or cell array of nx2 arrays) with x and y locations
%         of points to be transformed from img2 (reference) to img1.
% N                        : number of iterations, see 'doc imregdemons'. Default: 100
% AccumulatedFieldSmoothing: see 'doc imregdemons'. Default: 1.5
% PyramidLevels            : see 'doc imregdemons'. Default: 3
% DisplayWaitbar           : see 'doc imregdemons'. Default: false
% xvec     : vector with x-coordinates of image (2nd dimension of image).
%            Default: [1:imdim(2)]
% yvec     : vector with y-coordinates of image (1st dimension of image)
%            Default: [1:imdim(1)]
% method   : interpolation method for interpolating the displacement field.
%            Default: linear
% crop     : flag (true/false) determining whether to crop the image to the
%            locations where points are given to speed up registration.
%            Default: false
% padding  : if crop=true, padding determines how much of the image around
%            the points is included for cropping. Default 5 (same units as
%            xvec/yvec).
%
% OUTPUT
% 1st output : registered image
% 2nd output : transformed points
%
% Bart Bolsterlee, Neuroscience Research Australia
% March 2020

% Edits:
% 24/4, BB: removed option to transform label field.
% 24/4, BB: shifted coordinate system by half a pixel so that xvec/yvec
% indicate the center of a pixel (was edge of the pixel)
% 1/5, BB: changed interpolation so that values outside the image range are
% 0 (instead of NaN). This means that points outside the image range will
% be returned as the same value as the input (displacement = 0).

%% Parse inputs
p = inputParser;
addRequired(p,'img1')
addRequired(p,'img2')
% addRequired(p,'label1')
% addParameter(p,'label1',[])
addParameter(p,'points',[])

% Registration settings
addParameter(p,'N',100,@isscalar)
addParameter(p,'AccumulatedFieldSmoothing',1.5)
addParameter(p,'PyramidLevels',3)
addParameter(p,'DisplayWaitbar',false)
addParameter(p,'xvec',[])
addParameter(p,'yvec',[])
addParameter(p,'method','linear')
addParameter(p,'crop',false)
addParameter(p,'padding',5,@isscalar)
parse(p,img1,img2,varargin{:});

% Demons registration settings. See 'doc imregdemons' for details.
N       = p.Results.N; % number of iterations
afs     = p.Results.AccumulatedFieldSmoothing;
pl      = p.Results.PyramidLevels;
dw      = p.Results.DisplayWaitbar;

% Other settings.
crop    = p.Results.crop;
padding = p.Results.padding;
points  = p.Results.points; % points in fixed image coordinates
xvec    = p.Results.xvec;
yvec    = p.Results.yvec;
method  = p.Results.method; % interpolation method

%% Register images

% first input  : moving image
% second input : fixed image (reference)
% D            : displacement vectors from the fixed image grid to a
%                corresponding location in the moving image in image pixel
%                dimensions.
% img1_reg     : warped moving image (should match the fixed image)

% Build up vector of x and y location of pixel centres.
imdim = size(img2);
if isempty(xvec)
    xvec = 1 : imdim(2);    
end
if isempty(yvec)
    yvec = 1 : imdim(1);    
end
spacing = [xvec(2)-xvec(1) yvec(2)-yvec(1)];

if crop == true && ~isempty(points)
    % Register only regions in which points are provided.
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

% Call imregdemons to calculate the displacement field D that maps from the
% fixed to the moving image.
img1_reg = NaN(size(img2));
[D,I_reg] = imregdemons(img1(sel_y,sel_x),img2(sel_y,sel_x),N,...
    'AccumulatedFieldSmoothing',afs,...
    'PyramidLevels',pl,...
    'DisplayWaitbar',dw);
img1_reg(sel_y,sel_x,:) = I_reg;


% Transform the points
if isempty(points)
    points_tf = [];
else    
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
        % Calculate the new location of points given the displacement field
        % D.
        d = zeros(size(pts));
        d(:,1) = interp2(xvec(sel_x),yvec(sel_y),D(:,:,1),pts(:,1),pts(:,2),method,0); % x displacement in pixel coordinates
        d(:,2) = interp2(xvec(sel_x),yvec(sel_y),D(:,:,2),pts(:,1),pts(:,2),method,0); % y displacement in pixel coordinates

        % Update the new position vectors. The displacement field
        % predicted by imregdemons is in pixel coordinates, so it needs
        % to be converted back to image coordinates (flipped x/y axis,
        % and multiplied by pixel size).
        pnew = pts + d.*spacing;        
    end
%% Outputs
if nargout > 0
    % First output: registered image.
    varargout{1} = img1_reg;
    if nargout > 1
        % Second output: transformed points.
        varargout{2} = points_tf;
    end
end


end

