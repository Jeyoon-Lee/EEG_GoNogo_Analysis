% ERN 측정창 민감도 분석
%   여러 측정창에서 (1)조건효과 (2)집단x조건 상호작용이 일관된지 확인.
%   데이터(subject별 FA/Hit ERP)는 한 번만 로드하고, 측정창만 바꿔 amplitude 재계산.
%
% 전제: Compute_ERP.m 의 collect_group 과 동일한 방식으로 subject별 파형을 모은다.
%   (파형 = ROI 평균, baseline correction 적용된 [subject x time])

clear; close all; clc;

%% ===== PARAMETERS =====
proc_dir   = 'ERP/Data/T2/processed_data_resp_-400_600';
Meta_dir   = 'ERP/Raw_Data_Info/NESTdata_fromCCPL_260604.xlsx';
Pre_dir    = 'ERP/Data/T2/MADE_report_260704.csv';
roi_labels = {'E4','E7','E54'};
baseline_window = [-200 -100];
min_fa     = 4;

% 테스트할 측정창들 (ms)
windows = {[0 50], [-20 50], [0 60], [0 80], [0 100], [25 75], [40 80]};

%% ===== 대상 선정 =====
Meta_T = readtable(Meta_dir);
Pre_T  = readtable(Pre_dir);
Pre_T.ID = string(regexp(Pre_T.filename, 'NT\d{3}', 'match', 'once'));
[tf, loc] = ismember(Pre_T.ID, string(Meta_T.ID));
Pre_T.Diagnosis = nan(height(Pre_T),1);
Pre_T.Diagnosis(tf) = Meta_T.Diagnosis(loc(tf));
valid = Pre_T.n_epoch_nogo_fa >= min_fa & ~isnan(Pre_T.Diagnosis);
asd_files = Pre_T.filename(valid & Pre_T.Diagnosis==1);
td_files  = Pre_T.filename(valid & Pre_T.Diagnosis==0);
fprintf('대상: ASD %d, TD %d\n\n', numel(asd_files), numel(td_files));

%% ===== subject별 ERP 파형 수집 (한 번만) =====
[asd_err_wave, asd_cor_wave, times] = collect_waves(asd_files, proc_dir, roi_labels, baseline_window);
[td_err_wave,  td_cor_wave,  ~    ] = collect_waves(td_files,  proc_dir, roi_labels, baseline_window);

%% ===== 각 측정창에서 통계 =====
fprintf('%-14s %8s %8s %10s %8s %8s %8s\n', ...
    'Window(ms)', 'FA', 'Hit', 'Cond p', 'dz', 'Intx p', 'dGroup');
fprintf('%s\n', repmat('-',1,72));

results = [];
for k = 1:numel(windows)
    w = windows{k};
    widx = times >= w(1) & times <= w(2);

    % subject별 mean amplitude
    asd_fa  = mean(asd_err_wave(:, widx), 2);
    asd_hit = mean(asd_cor_wave(:, widx), 2);
    td_fa   = mean(td_err_wave(:,  widx), 2);
    td_hit  = mean(td_cor_wave(:,  widx), 2);

    all_fa  = [asd_fa;  td_fa];
    all_hit = [asd_hit; td_hit];

    % (Q1) 조건효과: paired t (FA vs Hit)
    [~, p_cond, ~, st_cond] = ttest(all_fa, all_hit);
    dz = mean(all_fa - all_hit) / std(all_fa - all_hit);   % paired Cohen's d

    % (Q3) 상호작용 = ΔERN 집단비교 (independent t)
    d_asd = asd_fa - asd_hit;
    d_td  = td_fa  - td_hit;
    [~, p_intx, ~, st_intx] = ttest2(d_asd, d_td);
    % 집단 효과크기 (Cohen's d, pooled)
    n1=numel(d_asd); n2=numel(d_td);
    sp = sqrt(((n1-1)*var(d_asd)+(n2-1)*var(d_td))/(n1+n2-2));
    d_group = (mean(d_asd)-mean(d_td))/sp;

    fprintf('[%4d %4d]    %7.3f %7.3f %10.2e %7.2f %8.4f %7.3f\n', ...
        w(1), w(2), mean(all_fa), mean(all_hit), p_cond, dz, p_intx, d_group);

    results(k).win = w; %#ok<SAGROW>
    results(k).p_cond = p_cond; results(k).dz = dz;
    results(k).p_intx = p_intx; results(k).d_group = d_group;
    results(k).asd_dERN = mean(d_asd); results(k).td_dERN = mean(d_td);
end

fprintf('%s\n', repmat('-',1,72));
fprintf('\n해석:\n');
fprintf(' - Cond p: 모든 창에서 매우 유의해야 (ERN 존재)\n');
fprintf(' - Intx p: 창마다 값이 다를 수 있음. 방향(dGroup 부호)이 일관되면 신뢰↑\n');
fprintf(' - dGroup 음수 = ASD ΔERN이 TD보다 덜 negative = ASD ERN 감소\n');

%% ===== 요약 그림: 측정창별 상호작용 =====
figure('Name','Window sensitivity','Position',[100 100 900 400]);

subplot(1,2,1); hold on;
win_labels = cellfun(@(w) sprintf('%d~%d', w(1),w(2)), windows, 'uni',0);
asd_d = [results.asd_dERN]; td_d = [results.td_dERN];
bar([asd_d; td_d]');
set(gca,'XTick',1:numel(windows),'XTickLabel',win_labels,'XTickLabelRotation',30);
ylabel('\DeltaERN (FA - Hit, \muV)'); legend({'ASD','TD'},'Location','best');
title('\DeltaERN by window'); grid on;

subplot(1,2,2); hold on;
p_intx_all = [results.p_intx];
bar(p_intx_all);
plot(get(gca,'xlim'), [0.05 0.05], 'r--');
set(gca,'XTick',1:numel(windows),'XTickLabel',win_labels,'XTickLabelRotation',30);
ylabel('Interaction p'); title('Interaction p by window'); grid on;
text(numel(windows)/2, 0.06, 'p=0.05', 'Color','r');


%% ===================== Local functions =====================
function [err_wave, cor_wave, times] = collect_waves(file_list, proc_dir, roi_labels, bl_win)
% subject별 조건 ERP(ROI 평균 파형)를 행으로 쌓아 반환 [subject x time]
    err_wave = []; cor_wave = []; times = [];
    for s = 1:numel(file_list)
        setfile = strrep(char(file_list{s}), '.mff', '_processed_data.set');
        try
            EEG = pop_loadset('filename', setfile, 'filepath', proc_dir);
        catch
            warning('로드 실패: %s', setfile); continue;
        end
        EEG = pop_rmbase(EEG, bl_win);
        if isempty(times), times = EEG.times; end

        [ie, ic] = get_cond_idx(EEG);
        if isempty(ie) || isempty(ic), continue; end

        roi = zeros(1,numel(roi_labels));
        for c = 1:numel(roi_labels)
            roi(c) = find(strcmpi({EEG.chanlocs.labels}, roi_labels{c}), 1);
        end
        err_wave = [err_wave; mean(mean(EEG.data(roi,:,ie),3),1)]; %#ok<AGROW>
        cor_wave = [cor_wave; mean(mean(EEG.data(roi,:,ic),3),1)]; %#ok<AGROW>
    end
end

function [ie, ic] = get_cond_idx(EEG)
    n = EEG.trials; gg = nan(1,n); ac = nan(1,n);
    for e = 1:n
        lat = cell2mat(EEG.epoch(e).eventlatency);
        i0 = find(lat==0,1); if isempty(i0), continue; end
        g = EEG.epoch(e).eventGoNogo;   if iscell(g), g = g{i0}; end
        a = EEG.epoch(e).eventAccuracy; if iscell(a), a = a{i0}; end
        gg(e)=g; ac(e)=a;
    end
    ie = find(gg==2 & ac==0);   % NoGo-FA
    ic = find(gg==1 & ac==1);   % Go-Hit
end
