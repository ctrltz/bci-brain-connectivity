function allplots_cortex_subplots(sa, data, colorlimits, cm, unit, smooth, varargin)
% ALLPLOTS_CORTEX_SUBPLOTS Source space visualizations with the New York
% Head model
%
% Parameters:
%   sa - source analysis structure
%   data - data to be plotted on the cortical surface (75K values)
%   colorlimits - data limits for the colormap
%   cm - colormap to be used
%   unit - label for the colormap
%   smooth - whether to use smoothed (if 1) or raw (if 0) coordinates
%   varargin - additional parameters are described below
%
% Varargin options:
% example: (...,'views', [0 5; 2 4], 'save', 1, 'namesave','mysources','saveformat','epsc')
%
%   views: a matrix with numbers representing different views to be
%       arranged in the subplot grid according to their position in the matrix
%
%       0 = ignore subplot
%       1 = left lateral, 2 = left hem medial
%       3 = rightlateral, 4 = right medial
%       5 = dorsal, 6 = dorsal horizontal
%       7 = ventral, 8 = ventral horizontal
%       9 = colorbar
%
%       negative values (-1, ..., -8) can be used to add colorbar to the
%       subplot
%
%   save : if = 1,  figures will be saved, default = 0
%   namesave : the name of the figures. default='source'
%   saveformat: format of the saved figure. without dot. default = 'epsc'
%   sources: plot sources as spheres
%   source_colors: provide color for all the sources to be plotted
%   cb_fontname: font to use for the colorbar
%   cb_fontsize: font size to use for the colorbar
%   cb_location: control the location of the colorbar
%   ax: provide axes where to put the plot (if not provided, axes are taken
%       automatically from a subplot grid based on the views)
%
% The function was authored by Stefan Haufe and modified by Mina Jamshidi 
% and Nikolai Kapralov. 

%% Default settings
set(0,'DefaultFigureColor',[1 1 1])
if smooth 
    vc = sa.cortex75K.vc_smooth;
    sm = '_smooth';
else
    vc = sa.cortex75K.vc;
    sm = '';
end
surface_pars = struct('alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, ...
  'showdirections', 0, 'colorbars', 0, 'dipnames', [], 'mymarkersize', 10, 'directions', [0 0 1 1 1 1], ...
  'printcbar', 1, 'userticks', []);

%% Handle varargin

%defaults
DoSave = 0;
saveformat = 'epsc';
savenameflag = 0;
saveName = 'source';
sources = {};
cb_fontsize = 10;
cb_fontname = [];
cb_location = 'eastoutside';
custom_ax = [];

% check varargin
if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for j = 1:2:(length(varargin)-1)
        if ~ischar (varargin{j})
            error(['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        switch lower (varargin{j})
            case 'source_colors'
                surface_pars.dotveccolors = varargin{j+1};
            case 'sources'
                sources = varargin{j+1};
            case 'views'
                views = varargin{j+1};
            case 'save'
                DoSave = varargin{j+1};
            case 'savename'
                savenameflag = 1;
                saveName = varargin{j+1};
            case 'saveformat'
                saveformat = varargin{j+1};
            case 'cb_fontsize'
                cb_fontsize = varargin{j+1};
            case 'cb_fontname'
                cb_fontname = varargin{j+1};
            case 'cb_location'
                cb_location = varargin{j+1};
            case 'ax'
                custom_ax = varargin{j+1};
            otherwise
                error(['Unknown name of the optional parameter: ' varargin{j}]);
        end
        
    end
end

if DoSave
    if ~savenameflag
        warning('No name for figures specified. The default name ((source)) is used')
    end
end

%% Plots
assert(ismatrix(views));
if ~isempty(custom_ax)
    assert(isscalar(views));
end

[rows, cols] = size(views);
% flatten the views matrix to iterate over subplots in the correct order
views_flat = views';
views_flat = views_flat(:);

for i = 1:numel(views_flat)
    % ignore 0s
    if ~views_flat(i)
        continue;
    end
    
    % negative values -> add colorbar
    make_colorbar = false;
    if views_flat(i) < 0
        views_flat(i) = -views_flat(i);
        make_colorbar = true;
    end

    % Use subplots if no axes were provided
    if isempty(custom_ax)
        ax = subplot(rows, cols, i);
    else
        ax = custom_ax;
    end
    
    switch (views_flat(i))
        case 1  %  view = 1, left lateral
            surface_pars.myviewdir = [-1 0 0];
            showsurface3(vc, sa.cortex75K.tri_left, surface_pars, data, sources{:});
        case 2  %  view = 2, left medial
            surface_pars.myviewdir = [-1 0 0];
            showsurface3(vc, sa.cortex75K.tri_right, surface_pars, data, sources{:});
        case 3  %  view = 3, right lateral
            surface_pars.myviewdir = [1 0 0];
            showsurface3(vc, sa.cortex75K.tri_right, surface_pars, data, sources{:});
        case 4  %  view = 4; right medial
            surface_pars.myviewdir = [1 0 0];
            showsurface3(vc, sa.cortex75K.tri_left, surface_pars, data, sources{:});
        case 5  %  view = 5, dorsal view
            surface_pars.myviewdir = [0 0 1];
            showsurface3(vc, sa.cortex75K.tri, surface_pars, data, sources{:});
        case 6  %  view = 6; dorsal view - horizontal
            surface_pars.myviewdir = [-1e-10 0 1];
            surface_pars.directions = [1 0 1 1 0 0];
            showsurface3(vc, sa.cortex75K.tri, surface_pars, data, sources{:}); %upperview rotated
        case 7  %  view = 7, ventral view
            surface_pars.myviewdir = [-1e-10 0 -1]; 
            showsurface3(vc, sa.cortex75K.tri, surface_pars, data, sources{:});
        case 8  %  view = 8; ventral view - horizontal
            surface_pars.myviewdir = [0 1e-10 -1];
            showsurface3(vc, sa.cortex75K.tri, surface_pars, data, sources{:});
        case 9  %  colorbar
            hf = imagesc(randn(5)); colormap(ax, cm)
            set(gca, 'clim', colorlimits, 'visible', 'off')
            set(hf, 'visible', 'off')
            cb = colorbar(cb_location); 
            set(cb, 'fontsize', cb_fontsize)
            if ~isempty(cb_fontname)
                set(cb, 'fontname', cb_fontname)
            end
            if ~isempty(surface_pars.userticks)
                set(cb, 'xtick', sort([colorlimits, surface_pars.userticks]))
            end
            ylabel(cb, unit)
        otherwise
            error(['Bad value in the views matrix ' num2str(views_flat(i))]);
    end

    if (make_colorbar)
        colormap(ax, cm);
        set(ax, 'clim', colorlimits, 'visible', 'off');
        cb = colorbar(cb_location); 
        set(cb, 'fontsize', cb_fontsize);
        if ~isempty(cb_fontname)
            set(cb, 'fontname', cb_fontname)
        end
        if ~isempty(surface_pars.userticks)
            set(cb, 'xtick', sort([colorlimits, surface_pars.userticks]));
        end
        ylabel(cb, unit);
    end
end

if DoSave
    saveas(gcf,saveName,saveformat);
end

set(0,'DefaultFigureColor',[1 1 1])

end
