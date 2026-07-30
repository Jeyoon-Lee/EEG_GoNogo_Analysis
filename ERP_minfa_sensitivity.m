% 스크리닝 기준 민감도 분석
%   측정창은 고정. 스크리닝 기준(min_fa / min_go / min_nogo)을 바꿔가며
%   표본 수와 통계 변화를 봄. 쓰지 않는 기준은 0으로 두면 무시
%     ERN/Pe (response-locked) : nogo_acc=0, 기준 [min_fa 0 0]
%     N2/P3  (stimulus-locked) : nogo_acc=1, 기준 [0 min_go min_nogo]
%   전체 subject 파형을 한 번만 로드하고, 필터만 바꿔 재계산

clear; close all; clc;

%% ===== PARAMETERS =====
proc_dir   = 'ERP/Data/T2/processed_data_resp_-400_600';
Meta_dir   = 'ERP/Raw_Data_Info/NESTdata_fromCCPL_260604.xlsx';
Pre_dir    = 'ERP/Data/T2/MMADE_resp_report_260724.csv'; % different dir for stm+/resp
roi_labels = {'E4','E7','E54'};
baseline_window = [-200 -100];       % 기준선 보정 구간
meas_window     = [-20 50];          % mean amplitude 측정 구간 (고정)

nogo_acc = 0;    % NoGo 조건 정의: 0 = NoGo-FA (ERN/Pe) | 1 = NoGo-CR (N2/P3)

% 테스트할 스크리닝 기준: 한 행이 [min_fa min_go min_nogo]
% N2/P3용
% crit_list = [ 0   20   10     
%               0  150   20   
%               0  150   30
%               0  150   40
%               0  150   50
%               0  150   60
%               0   20   40    
%               0  100   40
%               0  180   40
%               0  200   40
%               0  180   50 ];  
% ERN/Pe용
crit_list = [ 4  0  0
              6  0  0
              8  0  0
             10  0  0
             15  0  0 ];

%% ID - 진단 매핑
Meta_T = readtable(Meta_dir);
Pre_T  = readtable(Pre_dir);
Pre_T.ID = string(regexp(Pre_T.filename, 'NT\d{3}', 'match', 'once'));
[tf, loc] = ismember(Pre_T.ID, string(Meta_T.ID));
Pre_T.Diagnosis = nan(height(Pre_T),1);
Pre_T.Diagnosis(tf) = Meta_T.Diagnosis(loc(tf));

% 진단정보 있는 전체 대상
has_dx    = ~isnan(Pre_T.Diagnosis);
all_files = Pre_T.filename(has_dx);
all_dx    = Pre_T.Diagnosis(has_dx);
all_n_fa   = Pre_T.n_epoch_nogo_fa(has_dx);   % NoGo-FA  trial 수
all_n_go   = Pre_T.n_epoch_go_hit(has_dx);    % Go-Hit   trial 수
all_n_nogo = Pre_T.n_epoch_nogo_cr(has_dx);   % NoGo-CR  trial 수

%% 전체 subject 파형 수집
fprintf('전체 %d명 파형 로드 중...\n', numel(all_files));
[nogo_amp, go_amp, keep_dx, keep_n_fa, keep_n_go, keep_n_nogo] = deal([]);
widx_ready = false; times = [];
for i = 1:numel(all_files)
    setfile = strrep(char(all_files{i}), '.mff', '_processed_data.set');
    try
        EEG = pop_loadset('filename', setfile, 'filepath', proc_dir);
    catch
        warning('로드 실패: %s', setfile); continue;
    end
    EEG = pop_rmbase(EEG, baseline_window);
    if ~widx_ready
        times = EEG.times;
        widx  = times >= meas_window(1) & times <= meas_window(2);
        widx_ready = true;
    end

    [i_nogo, i_go] = get_cond_idx(EEG, nogo_acc);
    if isempty(i_nogo) || isempty(i_go), continue; end

    roi = zeros(1,numel(roi_labels));
    for c = 1:numel(roi_labels)
        roi(c) = find(strcmpi({EEG.chanlocs.labels}, roi_labels{c}), 1);
    end
    erp_nogo = mean(mean(EEG.data(roi,:,i_nogo),3),1);
    erp_go   = mean(mean(EEG.data(roi,:,i_go),  3),1);

    nogo_amp    = [nogo_amp;    mean(erp_nogo(widx))]; %#ok<AGROW>
    go_amp      = [go_amp;      mean(erp_go(widx))];   %#ok<AGROW>
    keep_dx     = [keep_dx;     all_dx(i)];            %#ok<AGROW>
    keep_n_fa   = [keep_n_fa;   all_n_fa(i)];          %#ok<AGROW>
    keep_n_go   = [keep_n_go;   all_n_go(i)];          %#ok<AGROW>
    keep_n_nogo = [keep_n_nogo; all_n_nogo(i)];        %#ok<AGROW>
    
end
fprintf('로드 완료: %d명\n\n', numel(nogo_amp));

%% 각 스크리닝 기준에서 통계
fprintf('측정창 고정: %d~%d ms | ROI %s | NoGo 조건: %s\n', ...
    meas_window(1), meas_window(2), strjoin(roi_labels,'+'), ...
    ternary(nogo_acc==0, 'NoGo-FA', 'NoGo-CR'));
fprintf('%-7s %-7s %-9s %6s %6s %6s %11s %8s %9s %8s\n', ...
    'min_fa','min_go','min_nogo','N','ASD','TD','Cond p','dz','Intx p','dGroup');
fprintf('%s\n', repmat('-',1,90));

for r = 1:size(crit_list,1)
    m_fa = crit_list(r,1); m_go = crit_list(r,2); m_nogo = crit_list(r,3);

    % 세 기준을 모두 만족 (0으로 둔 기준은 자동으로 무시됨)
    sel = keep_n_fa >= m_fa & keep_n_go >= m_go & keep_n_nogo >= m_nogo;

    dx   = keep_dx(sel);
    nogo = nogo_amp(sel);
    go   = go_amp(sel);

    n_asd = sum(dx==1); n_td = sum(dx==0);

    % 조건효과 (paired)
    [~, p_cond] = ttest(nogo, go);
    dz = mean(nogo-go)/std(nogo-go);

    % 상호작용 = Δ(NoGo-Go) 집단비교
    % nogo_acc = 0 -> NoGo-FA vs. Go-Hit
    % nogo_acc = 1 -> NoGO-CR vs. Go-Hit
    d_asd = nogo(dx==1) - go(dx==1);
    d_td  = nogo(dx==0) - go(dx==0);
    [~, p_intx] = ttest2(d_asd, d_td);
    n1=numel(d_asd); n2=numel(d_td);
    sp = sqrt(((n1-1)*var(d_asd)+(n2-1)*var(d_td))/(n1+n2-2));
    d_group = (mean(d_asd)-mean(d_td))/sp;

    fprintf('>=%-5d >=%-5d >=%-7d %6d %6d %6d %11.2e %8.2f %9.4f %8.3f\n', ...
        m_fa, m_go, m_nogo, numel(nogo), n_asd, n_td, p_cond, dz, p_intx, d_group);
end
fprintf('%s\n', repmat('-',1,90));
fprintf('\n해석:\n');
fprintf(' - 0으로 둔 기준은 필터링에 영향 없음\n');
fprintf(' - N/ASD/TD: 기준 높일수록 표본 감소 (검정력 하락)\n');
fprintf(' - Cond p: 모든 기준에서 유의해야 (성분 존재)\n');
fprintf(' - dGroup 부호가 기준 무관하게 일관되면 강건\n');
fprintf(' - Intx p: 표본(검정력)과 측정정밀도의 트레이드오프로 변동 가능\n');


%% ===== Helpers =====
function [i_nogo, i_go] = get_cond_idx(EEG, nogo_acc)
% nogo_acc = 0 -> NoGo-FA (오류)      : ERN / Pe
% nogo_acc = 1 -> NoGo-CR (억제 성공) : N2 / P3
    n = EEG.trials; gg = nan(1,n); ac = nan(1,n);
    for e = 1:n
        % eventlatency는 이벤트가 여러 개면 cell, 하나뿐이면 숫자 -> stimulus-locked에서 무응답 trial은 stm+ 하나뿐이라 숫자가 됨
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