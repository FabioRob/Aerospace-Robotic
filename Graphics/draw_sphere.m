function [fig] = draw_sphere(position, radius, col, filled)

% This function draws a circle at coordinate pos3D (vector) of radius 'r' and color 'col'
% Arguments:
% position - coordinate to draw the sphere (vector)
% radius - sphere radius
% col - color 
% filled - flag that will make the sphere filled or not


if filled == 1
	fig = scatter3(position(1), position(2), position(3), radius, col, 'filled');
else
	fig = scatter3(position(1), position(2), position(3), radius, col);
end

end