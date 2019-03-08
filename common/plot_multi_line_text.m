function h = plot_multi_line_text( cell_str, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
parser = inputParser;
addOptional(parser,'left',.2)
addOptional(parser,'bottom',.6)
addOptional(parser,'width',.2)
addOptional(parser,'height',.2)
addOptional(parser,'del',.22)

parse(parser, varargin{:});
left = parser.Results.left;
bottom = parser.Results.bottom;
width = parser.Results.width;
height = parser.Results.height;
del = parser.Results.del;

current_axes = gca;
set(current_axes,'fontsize',14)

h = axes('Position', [left bottom width height], 'Layer','top');
%bar(x,y), title('Bar Title')

options = {'Interpreter','latex','units','normalized','FontSize',14};

cell_str = flipud(cell_str(:));
K = length(cell_str);
base = .1; x_shift = .1;
for k = 1:K
    text(x_shift, base+(k-1)*del, cell_str{k}, options{:});
end

axis(h, 'off', 'tight')
%axes(current_axes)

end

