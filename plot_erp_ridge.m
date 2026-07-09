%% plot_erp_ridge.m
% 채널을 전(anterior)→후(posterior)로 정렬한 ridge plot
%   각 채널의 subject-평균 ERP(ASD+TD 합침)를 FA(error)/Hit(correct) 나란히.
%   목적: 특정 peak(예: parietal Pe)이 후두쪽에 국한되는지, 전-후로 어떻게 변하는지 확인.

clear; close all; clc;

%% ===== 설정 =====
proc_dir   = 'ERP/Data/T2/processed_data_resp_-400_600';
Meta_dir   = 'ERP/Raw_Data_Info/NESTdata_fromCCPL_260604.xlsx';
Pre_dir    = 'ERP/Data/T2/MADE_report_260704.csv';
baseline_window = [-200 -100];
lowpass_hz = 30;                 % ERP용 저역통과 (Compute_ERP와 동일)
min_fa     = 8;
sort_axis  = 'auto';             % 'X','Y', 또는 'auto' (전후축 자동 감지)

%% ===== 대상 (ASD+TD 합침) =====
Meta_T = readtable(Meta_dir);
Pre_T  = readtable(Pre_dir);
Pre_T.ID = string(regexp(Pre_T.filename, 'NT\d{3}', 'match', 'once'));
[tf, loc] = ismember(Pre_T.ID, string(Meta_T.ID));
Pre_T.Diagnosis = nan(height(Pre_T),1);
Pre_T.Diagnosis(tf) = Meta_T.Diagnosis(loc(tf));
valid = Pre_T.n_epoch_nogo_fa >= min_fa & ~isnan(Pre_T.Diagnosis);
files = Pre_T.filename(valid);
fprintf('대상: %d명 (ASD+TD 합침)\n', numel(files));

%% ===== subject별 조건 ERP 수집 (전채널) =====
% fa_all, cor_all: [subject x channel x time]
fa_all = []; cor_all = []; times = []; chanlocs = [];
for s = 1:numel(files)
    setfile = strrep(char(files{s}), '.mff', '_processed_data.set');
    try
        EEG = pop_loadset('filename', setfile, 'filepath', proc_dir);
    catch
        warning('로드 실패: %s', setfile); continue;
    end
    [~] = evalc(sprintf('EEG = pop_eegfiltnew(EEG, [], %d);', lowpass_hz));
    EEG = pop_rmbase(EEG, baseline_window);
    if isempty(times), times = EEG.times; chanlocs = EEG.chanlocs; end

    [ie, ic] = get_cond_idx(EEG);
    if isempty(ie) || isempty(ic), continue; end

    fa_all(end+1,:,:)  = mean(EEG.data(:,:,ie), 3); %#ok<SAGROW>
    cor_all(end+1,:,:) = mean(EEG.data(:,:,ic), 3); %#ok<SAGROW>
end
fprintf('수집 완료: %d명\n', size(fa_all,1));

% grand average (subject 평균) → [channel x time]
ga_fa  = squeeze(mean(fa_all, 1));
ga_cor = squeeze(mean(cor_all, 1));

%% ===== 채널 전후축 정렬 =====
nb = numel(chanlocs);
X = [chanlocs.X];  Y = [chanlocs.Y];
% 전후축 자동 감지: 값의 분산이 큰 축을 앞뒤로 (EEGLAB에서 X=코방향이 보통 전후)
if strcmpi(sort_axis,'auto')
    if range(X) >= range(Y), aval = X; usedax = 'X';
    else,                    aval = Y; usedax = 'Y';
    end
elseif strcmpi(sort_axis,'X'), aval = X; usedax = 'X';
else,                          aval = Y; usedax = 'Y';
end
% 앞(anterior, 큰 X)이 위로 오도록 내림차순 정렬
[~, order] = sort(aval, 'descend');
fprintf('정렬축: %s (앞→뒤)\n', usedax);

%% ===== Ridge plot =====
figure('Name','ERP ridge (anterior→posterior)','Position',[100 60 800 1100]);
hold on;

% 각 채널을 세로로 offset해서 겹치지 않게
gap = 4;   % 채널 간 세로 간격 (µV 단위, 데이터 크기에 맞게 조정)
n_show = nb;
yticks_pos = zeros(1, n_show);
yticks_lab = cell(1, n_show);

for r = 1:n_show
    ch = order(r);
    yoff = (n_show - r) * gap;    % 위쪽이 anterior
    yticks_pos(r) = yoff;
    yticks_lab{r} = chanlocs(ch).labels;

    plot(times, ga_fa(ch,:)  + yoff, 'r', 'LineWidth', 0.8);   % FA (error)
    % plot(times, ga_cor(ch,:) + yoff, 'b', 'LineWidth', 0.8);   % Hit (correct)
    plot([times(1) times(end)], [yoff yoff], 'Color',[.8 .8 .8], 'LineWidth',0.3);  % 채널 기준선
end

plot([0 0], [min(yticks_pos)-gap, max(yticks_pos)+gap], 'k:', 'LineWidth', 1);  % 반응 시점
set(gca, 'YTick', fliplr(yticks_pos), 'YTickLabel', fliplr(yticks_lab), 'FontSize', 7);
xlabel('Time (ms)'); ylabel('Channel (anterior → posterior, top→bottom)');
title(sprintf('Grand-average ERP ridge  |  red=FA(error), blue=Hit(correct)  |  n=%d', size(fa_all,1)));
xlim([times(1) times(end)]);
ylim([min(yticks_pos)-gap, max(yticks_pos)+gap]);

legend_h(1) = plot(nan,nan,'r','LineWidth',1.5);
% legend_h(2) = plot(nan,nan,'b','LineWidth',1.5);
% legend(legend_h, {'FA (error)','Hit (correct)'}, 'Location','northeast');

fprintf('\ngap=%.1f µV로 채널 간격 설정. 파형이 겹치거나 너무 납작하면 gap 조정.\n', gap);


%% ===== Local functions =====
function [ie, ic] = get_cond_idx(EEG)
    n = EEG.trials; gg = nan(1,n); ac = nan(1,n);
    for e = 1:n
        lat = cell2mat(EEG.epoch(e).eventlatency);
        i0 = find(lat==0,1); if isempty(i0), continue; end
        g = EEG.epoch(e).eventGoNogo;   if iscell(g), g = g{i0}; end
        a = EEG.epoch(e).eventAccuracy; if iscell(a), a = a{i0}; end
        gg(e)=g; ac(e)=a;
    end
    ie = find(gg==2 & ac==0);
    ic = find(gg==1 & ac==1);
end
