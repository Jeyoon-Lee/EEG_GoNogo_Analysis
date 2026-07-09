function [cfg] = plot_raw(EEG, highlight_chans)
% PLOT_RAW visual inspection of EEGLAB data using fieldtrip ft_databrowser
% function
    if nargin <2
        highlight_chans = {};
    end

    data = eeglab2fieldtrip(EEG, 'raw', 'none');
    event = struct('type', {EEG.event.type}, ...
                   'sample', num2cell([EEG.event.latency]), ...
                   'value', {EEG.event.type}, ...
                   'offset', 0, ...
                   'duration', 0);
    
    cfg = [];
    cfg.continuous = 'yes';
    cfg.viewmode   = 'vertical';
    cfg.blocksize  = 5;
    cfg.event      = event;    
    cfg.plotlabels = 'yes'; 
    cfg.fontsize   = 8;
    cfg.ylim       = [-30 30];
       
    if ~isempty(highlight_chans)
        grp = ones(numel(data.label), 1);
        grp(ismember(data.label, highlight_chans)) = 2;
        cfg.colorgroups = grp;
        cfg.channelcolormap = [0 0 0; 1 0 0]; % group 1 = black; group 2 = red
    end
    
    [cfg] = ft_databrowser(cfg, data); 
end

