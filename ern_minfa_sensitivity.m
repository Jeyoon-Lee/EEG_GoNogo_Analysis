%% ern_minfa_sensitivity.m
% min_fa(최소 FA trial 기준) 민감도 분석
%   측정창은 -20~50 ms 고정. min_fa를 바꿔가며 표본 수와 통계 변화를 본다.
%   전체 subject 파형을 한 번만 로드하고, min_fa 필터만 바꿔 재계산.

clear; close all; clc;

%% ===== 설정 =====
proc_dir   = 'ERP/Data/T2/processed_data_resp_-400_600';
Meta_dir   = 'ERP/Raw_Data_Info/NESTdata_fromCCPL_260604.xlsx';
Pre_dir    = 'ERP/Data/T2/MADE_report_260704.csv';
roi_labels = {'E4','E7','E54'};
baseline_window = [-200 -100];
ern_window = [-20 50];               % 고정
minfa_list = [4 6 8 10 15];          % 테스트할 기준

%% ===== 메타 매칭 (전체) =====
Meta_T = readtable(Meta_dir);
Pre_T  = readtable(Pre_dir);
Pre_T.ID = string(regexp(Pre_T.filename, 'NT\d{3}', 'match', 'once'));
[tf, loc] = ismember(Pre_T.ID, string(Meta_T.ID));
Pre_T.Diagnosis = nan(height(Pre_T),1);
Pre_T.Diagnosis(tf) = Meta_T.Diagnosis(loc(tf));

% 진단정보 있는 전체 대상 (FA 기준은 나중에 필터)
has_dx = ~isnan(Pre_T.Diagnosis);
all_files = Pre_T.filename(has_dx);
all_dx    = Pre_T.Diagnosis(has_dx);
all_fa    = Pre_T.n_epoch_nogo_fa(has_dx);

%% ===== 전체 subject 파형 수집 (한 번만) =====
% 각 subject의 FA/Hit ERP(ROI 평균, -20~50 mean amp까지) + FA수 + 진단 저장
fprintf('전체 %d명 파형 로드 중...\n', numel(all_files));
[fa_amp, hit_amp, keep_dx, keep_fa] = deal([]);
widx_ready = false; times = [];
for i = 1:numel(all_files)
    setfile = strrep(char(all_files{i}), '.mff', '_processed_data.set');
    try
        EEG = pop_loadset('filename', setfile, 'filepath', proc_dir);
    catch
        warning('로드 실패: %s', setfile); continue;
    end
    EEG = pop_rmbase(EEG, baseline_window);
    if ~widx_ready, times = EEG.times; widx = times>=ern_window(1) & times<=ern_window(2); widx_ready = true; end

    [ie, ic] = get_cond_idx(EEG);
    if isempty(ie) || isempty(ic), continue; end

    roi = zeros(1,numel(roi_labels));
    for c = 1:numel(roi_labels)
        roi(c) = find(strcmpi({EEG.chanlocs.labels}, roi_labels{c}), 1);
    end
    erp_fa  = mean(mean(EEG.data(roi,:,ie),3),1);
    erp_hit = mean(mean(EEG.data(roi,:,ic),3),1);

    fa_amp  = [fa_amp;  mean(erp_fa(widx))];   %#ok<AGROW>
    hit_amp = [hit_amp; mean(erp_hit(widx))];  %#ok<AGROW>
    keep_dx = [keep_dx; all_dx(i)];            %#ok<AGROW>
    keep_fa = [keep_fa; all_fa(i)];            %#ok<AGROW>
end
fprintf('로드 완료: %d명\n\n', numel(fa_amp));

%% ===== 각 min_fa 기준에서 통계 =====
fprintf('측정창 고정: %d~%d ms | ROI %s\n', ern_window(1),ern_window(2), strjoin(roi_labels,'+'));
fprintf('%-8s %6s %6s %6s %11s %8s %9s %8s\n', ...
    'min_fa','N','ASD','TD','Cond p','dz','Intx p','dGroup');
fprintf('%s\n', repmat('-',1,66));

for m = minfa_list
    sel = keep_fa >= m;
    dx  = keep_dx(sel);
    fa  = fa_amp(sel);
    hit = hit_amp(sel);

    n_asd = sum(dx==1); n_td = sum(dx==0);

    % 조건효과 (paired)
    [~, p_cond] = ttest(fa, hit);
    dz = mean(fa-hit)/std(fa-hit);

    % 상호작용 = ΔERN 집단비교
    d_asd = fa(dx==1) - hit(dx==1);
    d_td  = fa(dx==0) - hit(dx==0);
    [~, p_intx] = ttest2(d_asd, d_td);
    n1=numel(d_asd); n2=numel(d_td);
    sp = sqrt(((n1-1)*var(d_asd)+(n2-1)*var(d_td))/(n1+n2-2));
    d_group = (mean(d_asd)-mean(d_td))/sp;

    fprintf('>=%-6d %6d %6d %6d %11.2e %8.2f %9.4f %8.3f\n', ...
        m, numel(fa), n_asd, n_td, p_cond, dz, p_intx, d_group);
end
fprintf('%s\n', repmat('-',1,66));
fprintf('\n해석:\n');
fprintf(' - N/ASD/TD: 기준 높일수록 표본 감소 (검정력 하락)\n');
fprintf(' - Cond p: 모든 기준에서 유의해야 (ERN 존재)\n');
fprintf(' - dGroup 양수 = ASD ERN 감소. 기준 무관하게 방향 일관되면 강건.\n');
fprintf(' - Intx p: 표본(검정력)과 측정정밀도의 트레이드오프로 변동 가능.\n');


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
