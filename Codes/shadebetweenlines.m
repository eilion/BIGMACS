function shadebetweenlines(x, y_upper,y_lower, cc)

% Shade between two lines with specified color

for qq = 2:length(y_upper)
    y1 = y_upper(qq-1);
    y2 = y_lower(qq-1);
    y3 = y_lower(qq);
    y4 = y_upper(qq);
    x1 = x(qq-1);
    x2 = x(qq-1);
    x3 = x(qq);
    x4 = x(qq);
    
    patch([x1 x2 x3 x4], [y1 y2 y3 y4], cc, 'FaceAlpha', 0.4, 'EdgeColor', 'none','handlevisibility','off')
    
end
