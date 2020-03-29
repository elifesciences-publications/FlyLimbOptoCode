function ConfAxis(fontSize)
if nargin < 1, fontSize = 12; end
set(gca, 'FontSize', fontSize, 'LineWidth',2, 'Box', 'Off');
end