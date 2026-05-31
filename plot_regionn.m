%% for countour overlap (misfit plane) plotting
function [ax,img]=plot_regionn(X, Y, data,min,max,color,lineStyle)
   
    pos = [0.1, 0.1, 0.8, 0.8]; 
    fig = figure('Visible', 'off');
    ax = axes('Position', pos, 'Parent', fig);
    [C,h]=contourf(ax,X,Y,data,[min,max]);

     % =
    cmap = [color;[1 1 1]];  % 
    colormap(gca, cmap);
    hold on;

    % 
    contour(X, Y, data, [min min],...
    'LineColor', color,...
    'LineStyle', lineStyle,...
    'LineWidth', 1.5);

    contour(X, Y, data, [max max],...
    'LineColor', color,...
    'LineStyle', lineStyle,...
    'LineWidth', 1.5);

    hFills = h.FacePrims; 
    set(hFills,'FaceAlpha',1);
    grid on;

    F1 = getframe(ax);
    img = F1.cdata;
    close(fig);

%     contour(X, Y, mask, [1 1],...
%         'LineColor', color,...
%         'LineStyle', lineStyle,...
%         'LineWidth', 1.5);
end