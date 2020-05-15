function PlotData = getPlotData3D(ObjectData,nLabels,zvec)
%%GETPLOTDATA3D creates a cell array with all 3D coordinates per label for 
% plotting purposes.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% April 2020

% Make plot coordinates per label
% nLabels = size(colors,1);
PlotData = cell(nLabels,1);
for slice_nr = 1 : size(ObjectData,1)
    for obj_nr = 1 : length(ObjectData(slice_nr).position)
        pos = ObjectData(slice_nr).position{obj_nr};
        type = ObjectData(slice_nr).type{obj_nr};
        label_nr = ObjectData(slice_nr).label_nr(obj_nr);
        switch type
            case 'spline'
                pos = fit_spline(pos);
            case 'polygon'
                pos = [pos;pos(1,:)];
        end
        % Add z-coordinate.
        pos(:,3) = zvec(slice_nr)*ones(size(pos,1),1);
        PlotData{label_nr} = [PlotData{label_nr};NaN(1,3);pos];
    end
end

end % of function