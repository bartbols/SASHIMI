function data = convertObjects(old_filename,img_filename,new_filename)
%CONVERTOBJECTS This script will convert object files created in version
%1.0 of SASHIMI Segmentation to later versions in which the x/y axis are
%flipped and shifted by half a pixel.

if nargin == 2
    new_filename = [];
end
% Load the object and the NIfTI header.
if nargin == 0
    % If no inputs are given, select files manually.
    % Old object file
    [file,path] = uigetfile('*.mat','Select the old object file');
    if file==0
        return
    else
        old_filename = fullfile(path,file);
    end
    
    % Corresponding image file.
    [file,path] = uigetfile({'*.nii;*.nii.gz;',...
                        'Image files (*.nii,*.nii.gz)'},...
                        'Select the corresponding image');
    if file==0
        return
    else
        img_filename = fullfile(path,file);
    end
    
    % New object file. 
    defname = [old_filename(1:end-4) '_converted.mat'];
    [file,path] = uiputfile('*.mat','Set filename of the converted object file',defname);
    if file==0
        return
    else
        new_filename = fullfile(path,file);
    end
end
data     = load(old_filename);
metadata = niftiinfo(img_filename);

%% Convert the object data to the new coordinate system used in SASHIMI 
% v1.1 and up.
OldObjectData    = data.ObjectData;
NewObjectData    = OldObjectData;

for slice_nr = 1 : length(OldObjectData)
%     NewObjectData(slice_nr).label_nr = OldObjectData(slice_nr).label_nr;
%     NewObjectData(slice_nr).type = OldObjectData(slice_nr).type;
    for i = 1 : length(OldObjectData(slice_nr).position)
        pos = OldObjectData(slice_nr).position{i};
        NewObjectData(slice_nr).position{i} = pos(:,[2 1]) - metadata.PixelDimensions(1:2)/2;
    end
    
end
data.ObjectData = NewObjectData;

%% Save the objects to a new file.
if ~isempty(new_filename)
    save(new_filename,'-struct','data');
    fprintf('Objects were successfully converted.\nOld object file: %s\nNew object file: %s\n',old_filename,new_filename)
end

end % of function

