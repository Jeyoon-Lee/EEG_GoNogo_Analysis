% ERP 측정창 민감도 분석
%   여러 측정창에서 (1)조건효과 (2)집단x조건 상호작용이 일관된지 확인
%   데이터(subject별 파형)는 한 번만 로드하고, 측정창만 바꿔 amplitude 재계산
%
%   스크리닝 기준(min_fa / min_go / min_nogo)은 쓰지 않는 것을 0으로 두면 무시
%     ERN/Pe (response-locked) : nogo_acc=0, min_fa=8, min_go=0,  min_nogo=0
%     N2/P3  (stimulus-locked) : nogo_acc=1, min_fa=0, min_go=40, min_nogo=15

clear; close all; clc;

%% ===== PARAMETERS =====
proc_dir   = 'ERP/Data/T2/processed_data_resp_-400_600';
Meta_dir   = 'ERP/Raw_Data_Info/NESTdata_fromCCPL_260604.xlsx';
Pre_dir    = 'ERP/Data/T2/MADE_report_260704.csv';
roi_labels = {'E4','E7','E54'};
baseline_window = [-200 -100];       % 기준선 보정 구간

nogo_acc = 0;    % NoGo 조건 정의: 0 = NoGo-FA (ERN/Pe) | 1 = NoGo-CR (N2/P3)

% 스크리닝 기준 (쓰지 않는 것은 0)
min_fa   = 8;
min_go   = 0;
min_nogo = 0;

% 테스트할 측정창들 (ms)
windows = {[0 50], [-20 50], [0 60], [0 80], [0 100], [25 75], [40 80]};

%% 스크리닝 기준에 맞는 파일 로드
Meta_T = readtable(Meta_dir);
Pre_T  = readtable(Pre_dir);
Pre_T.ID = string(regexp(Pre_T.filename, 'NT\d{3}', 'match', 'once'));
[tf, loc] = ismember(Pre_T.ID, string(Meta_T.ID));
Pre_T.Diagnosis = nan(height(Pre_T),1);
Pre_T.Diagnosis(tf) = Meta_T.Diagnosis(loc(tf));

% 세 기준을 모두 만족 (0으로 둔 기준은 자동으로 무시됨)
valid = Pre_T.n_epoch_nogo_fa >= min_fa ...
      & Pre_T.n_epoch_go_hit  >= min_go ...
      & Pre_T.n_epoch_nogo_cr >= min_nogo ...
      & ~isnan(Pre_T.Diagnosis);

asd_files = Pre_T.filename(valid & Pre_T.Diagnosis==1);
td_files  = Pre_T.filename(valid & Pre_T.Diagnosis==0);
fprintf('스크리닝 min_fa>=%d, min_go>=%d, min_nogo>=%d | NoGo 조건: %s\n', ...
    min_fa, min_go, min_nogo, ternary(nogo_acc==0,'NoGo-FA','NoGo-CR'));
fprintf('대상: ASD %d, TD %d\n\n', numel(asd_files), numel(td_files));

% subject별 ERP 파형 수집
[asd_nogo_wave, asd_go_wave, times] = collect_waves(asd_files, proc_dir, roi_labels, baseline_window, nogo_acc);
[td_nogo_wave,  td_go_wave,  ~    ] = collect_waves(td_files,  proc_dir, roi_labels, baseline_window, nogo_acc);

%% 각 측정창에서 통계
fprintf('%-14s %8s %8s %10s %8s %8s %8s\n', ...
    'Window(ms)', 'NoGo', 'Go', 'Cond p', 'dz', 'Intx p', 'dGroup');
fprintf('%s\n', repmat('-',1,72));

results = [];
for k = 1:numel(windows)
    w = windows{k};
    widx = times >= w(1) & times <= w(2);

    % subject별 mean amplitude
    asd_nogo = mean(asd_nogo_wave(:, widx), 2);
    asd_go   = mean(asd_go_wave(:,   widx), 2);
    td_nogo  = mean(td_nogo_wave(:,  widx), 2);
    td_go    = mean(td_go_wave(:,    widx), 2);

    all_nogo = [asd_nogo; td_nogo];
    all_go   = [asd_go;   td_go];

    % (Q1) 조건효과: paired t
    [~, p_cond] = ttest(all_nogo, all_go);
    dz = mean(all_nogo - all_go) / std(all_nogo - all_go);   % paired Cohen's d

    % (Q3) 상호작용 = Δ(NoGo-Go) 집단비교 (independent t)
    d_asd = asd_nogo - asd_go;
    d_td  = td_nogo  - td_go;
    [~, p_intx] = ttest2(d_asd, d_td);
    n1=numel(d_asd); n2=numel(d_td);
    sp = sqrt(((n1-1)*var(d_asd)+(n2-1)*var(d_td))/(n1+n2-2));
    d_group = (mean(d_asd)-mean(d_td))/sp;

    fprintf('[%4d %4d]    %7.3f %7.3f %10.2e %7.2f %8.4f %7.3f\n', ...
        w(1), w(2), mean(all_nogo), mean(all_go), p_cond, dz, p_intx, d_group);

    results(k).win = w; %#ok<SAGROW>
    results(k).p_cond = p_cond; results(k).dz = dz;
    results(k).p_intx = p_intx; results(k).d_group = d_group;
    results(k).asd_d  = mean(d_asd); results(k).td_d = mean(d_td);
end

fprintf('%s\n', repmat('-',1,72));
fprintf('\n해석:\n');
fprintf(' - Cond p: 모든 창에서 매우 유의해야 (성분 존재)\n');
fprintf(' - Intx p: 창마다 값이 다를 수 있음. 방향(dGroup 부호)이 일관되면 신뢰↑\n');
fprintf(' - dGroup 양수 = ASD Δ가 TD보다 덜 negative = ASD 성분 감소\n');

%% 요약 그림: 측정창별 상호작용
figure('Name','Window sensitivity','Position',[100 100 900 400]);

subplot(1,2,1); hold on;
win_labels = cellfun(@(w) sprintf('%d~%d', w(1),w(2)), windows, 'uni',0);
bar([[results.asd_d]; [results.td_d]]');
set(gca,'XTick',1:numel(windows),'XTickLabel',win_labels,'XTickLabelRotation',30);
ylabel('\Delta (NoGo - Go, \muV)'); legend({'ASD','TD'},'Location','best');
title('\Delta by window'); grid on;

subplot(1,2,2); hold on;
bar([results.p_intx]);
plot(get(gca,'xlim'), [0.05 0.05], 'r--');
set(gca,'XTick',1:numel(windows),'XTickLabel',win_labels,'XTickLabelRotation',30);
ylabel('Interaction p'); title('Interaction p by window'); grid on;
text(numel(windows)/2, 0.06, 'p=0.05', 'Color','r');


%% ===== Helpers =====
function [nogo_wave, go_wave, times] = collect_waves(file_list, proc_dir, roi_labels, bl_win, nogo_acc)
% subject별 조건 ERP(ROI 평균 파형)를 행으로 쌓아 반환 [subject x time]
    nogo_wave = []; go_wave = []; times = [];
    for s = 1:numel(file_list)
        setfile = strrep(char(file_list{s}), '.mff', '_processed_data.set');
        try
            EEG = pop_loadset('filename', setfile, 'filepath', proc_dir);
        catch
            warning('로드 실패: %s', setfile); continue;
        end
        EEG = pop_rmbase(EEG, bl_win);
        if isempty(times), times = EEG.times; end

        [i_nogo, i_go] = get_cond_idx(EEG, nogo_acc);
        if isempty(i_nogo) || isempty(i_go), continue; end

        roi = zeros(1,numel(roi_labels));
        for c = 1:numel(roi_labels)
            roi(c) = find(strcmpi({EEG.chanlocs.labels}, roi_labels{c}), 1);
        end
        nogo_wave = [nogo_wave; mean(mean(EEG.data(roi,:,i_nogo),3),1)]; %#ok<AGROW>
        go_wave   = [go_wave;   mean(mean(EEG.data(roi,:,i_go),  3),1)]; %#ok<AGROW>
    end
end

function [i_nogo, i_go] = get_cond_idx(EEG, nogo_acc)
% nogo_acc = 0 -> NoGo-FA (오류)      : ERN / Pe
% nogo_acc = 1 -> NoGo-CR (억제 성공) : N2 / P3
    n = EEG.trials; gg = nan(1,n); ac = nan(1,n);
    for e = 1:n
        % eventlatency는 이벤트가 여러 개면 cell, 하나뿐이면 숫자
        lat_raw = EEG.epoch(e).eventlatency;
        if iscell(lat_raw), lat = cell2mat(lat_raw); else, lat = lat_raw; end

        i0 = find(lat==0,1); if isempty(i0), continue; end
        g = EEG.epoch(e).eventGoNogo;   if iscell(g), g = g{i0}; end
        a = EEG.epoch(e).eventAccuracy; if iscell(a), a = a{i0}; end
        gg(e)=g; ac(e)=a;
    end
    i_nogo = find(gg==2 & ac==nogo_acc);   % NoGo (FA 또는 CR)
    i_go   = find(gg==1 & ac==1);          % Go-Hit
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end