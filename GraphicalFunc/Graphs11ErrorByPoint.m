function [figs] = Graphs11ErrorByPoint(pointTestForError, errorMatMultiscale,...
                            errorMatSinglescale, nLevels, nDim)
   figs = {};
   if nDim == 1
       [figs] = PlotGraph11_1D(pointTestForError, errorMatMultiscale,...
           errorMatSinglescale, nLevels);
   end
       
     
            
end

function [figs] = PlotGraph11_1D(pointTestForError, errorMatMultiscale,...
        errorMatSinglescale, nLevels)
figs = cell(nLevels, 1);
for i = 1 : size(errorMatMultiscale , 3)
   fig = figure;
   leg = {};
   colors = zeros(3, nLevels);
   for iLevel = 1 : nLevels
       hold on;
       ff  = plot(pointTestForError, errorMatMultiscale(:, iLevel, i), '-');
       leg{iLevel} = num2str(iLevel);
       colors(:, iLevel) = ff.Color;
   end
   for iLevel = 1 : nLevels
       plot(pointTestForError, errorMatSinglescale(:, iLevel, i), '--', 'Color', colors(:, iLevel));
   end
    legend(leg)
   switch i
       case 1
           title('Mean Error');
       case 2
           title('Median Error');
       case 3
           title('Max Error');
   end
   figs{i} = fig;
end



end