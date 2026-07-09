%% plot_ern_analysis.m
% ERN grand average 파형 + mean amplitude 추출 + 통계
%   채널: E4,E7,E54 평균 (frontocentral ROI)
%   조건: NoGo-FA(error) vs Go-Hit(correct)
%   ERN 정량화: mean amplitude, 0~100ms
%   baseline: -200~-100ms (Wang 2020)
%   포함: n_epoch_nogo_fa >= 6
%
% 통계:
%   Q1. 조건 효과 (FA vs Hit, 집단무관)         -> paired t-test (양측)
%   Q2. ERN < CRN 방향성 (FA가 더 negative)     -> paired t-test (단측)
%   Q3. 집단(ASD/TD) x 조건(FA/Hit) 상호작용    -> mixed ANOVA (fitrm/ranova)

clear; close all; clc;

%% ===== 설정 =====
proc_dir   = 'ERP/Data/T2/processed_data_resp_-400_600';
Meta_dir   = 'ERP/Raw_Data_Info/NESTdata_fromCCPL_260604.xlsx';
Pre_dir    = 'ERP/Data/T2/MADE_report_260704.csv';
roi_labels = {'E33', 'E34', 'E36', 'E38'}; % for ERN {'E4','E7','E54'}; for Pe {'E33', 'E34', 'E36', 'E38'};
baseline_window = [-200 -100];   
ern_window      = [200 400];       % ms, mean amplitude 측정 구간 % for ERN [-20 50]; for Pe [200 400]
min_fa     = 8;

%% ===== 메타 매칭 + 포함 대상 =====
Meta_T = readtable(Meta_dir);
Pre_T  = readtable(Pre_dir);
Pre_T.ID = string(regexp(Pre_T.filename, 'NT\d{3}', 'match', 'once'));
Meta_ID  = string(Meta_T.ID);
[tf, loc] = ismember(Pre_T.ID, Meta_ID);
Pre_T.Diagnosis = nan(height(Pre_T), 1);
Pre_T.Diagnosis(tf) = Meta_T.Diagnosis(loc(tf));

valid = Pre_T.n_epoch_nogo_fa >= min_fa & ~isnan(Pre_T.Diagnosis);
asd_files = Pre_T.filename(valid & Pre_T.Diagnosis == 1);
td_files  = Pre_T.filename(valid & Pre_T.Diagnosis == 0);
fprintf('포함: ASD %d명, TD %d명 (FA>=%d)\n', numel(asd_files), numel(td_files), min_fa);

%% ===== subject별 ERP + mean amplitude 수집 =====
[asd_err, asd_cor, asd_amp_err, asd_amp_cor, times] = ...
    collect_group(asd_files, proc_dir, roi_labels, baseline_window, ern_window);
[td_err,  td_cor,  td_amp_err,  td_amp_cor,  ~] = ...
    collect_group(td_files,  proc_dir, roi_labels, baseline_window, ern_window);

% grand average + SEM
ga = @(M) mean(M,1);
se = @(M) std(M,0,1)/sqrt(size(M,1));

%% ===== 파형 그리기 =====
figure('Name','ERN grand average','Position',[80 80 1400 450]);

subplot(1,3,1); hold on;
h1 = shaded(times, ga(asd_err), se(asd_err), [0.8 0.2 0.2]);
h2 = shaded(times, ga(asd_cor), se(asd_cor), [0.2 0.2 0.8]);
decorate(sprintf('ASD (n=%d)', numel(asd_files)), ern_window);
legend([h1 h2], {'Error(FA)','Correct(Hit)'}, 'Location','best');

subplot(1,3,2); hold on;
h1 = shaded(times, ga(td_err), se(td_err), [0.8 0.2 0.2]);
h2 = shaded(times, ga(td_cor), se(td_cor), [0.2 0.2 0.8]);
decorate(sprintf('TD (n=%d)', numel(td_files)), ern_window);
legend([h1 h2], {'Error(FA)','Correct(Hit)'}, 'Location','best');

subplot(1,3,3); hold on;
h1 = shaded(times, ga(asd_err-asd_cor), se(asd_err-asd_cor), [0.85 0.33 0.10]);
h2 = shaded(times, ga(td_err-td_cor),  se(td_err-td_cor),  [0.00 0.45 0.74]);
decorate('Difference (Error - Correct)', ern_window);
legend([h1 h2], {'ASD','TD'}, 'Location','best');

sgtitle(sprintf('ERN @ %s | baseline %d~%d | mean amp %d~%d ms', ...
    strjoin(roi_labels,'+'), baseline_window(1), baseline_window(2), ...
    ern_window(1), ern_window(2)));

%% ===== 통계 =====
% 전체(집단 무관) amplitude 벡터
amp_err = [asd_amp_err; td_amp_err];   % FA (error) mean amp, subject별
amp_cor = [asd_amp_cor; td_amp_cor];   % Hit (correct)
group   = [ones(numel(asd_amp_err),1); zeros(numel(td_amp_err),1)];  % 1=ASD,0=TD

fprintf('\n========== 통계 ==========\n');
fprintf('mean amplitude 측정: %d~%d ms, ROI %s\n\n', ...
    ern_window(1), ern_window(2), strjoin(roi_labels,'+'));

% --- Q1. 조건 효과 (FA vs Hit), 양측 paired t ---
[~, p1, ci1, st1] = ttest(amp_err, amp_cor);   % 양측
fprintf('[Q1] FA vs Hit (집단무관, paired t, 양측)\n');
fprintf('     FA  = %.3f ± %.3f µV\n', mean(amp_err), std(amp_err)/sqrt(numel(amp_err)));
fprintf('     Hit = %.3f ± %.3f µV\n', mean(amp_cor), std(amp_cor)/sqrt(numel(amp_cor)));
fprintf('     t(%d)=%.3f, p=%.4g\n\n', st1.df, st1.tstat, p1);

% --- Q2. ERN < CRN (FA가 더 negative), 단측 paired t ---
[~, p2, ~, st2] = ttest(amp_err, amp_cor, 'Tail','left');   % FA < Hit
fprintf('[Q2] ERN < CRN 방향성 (FA < Hit, paired t, 단측)\n');
fprintf('     t(%d)=%.3f, p=%.4g\n', st2.df, st2.tstat, p2);
d_paired = mean(amp_err-amp_cor)/std(amp_err-amp_cor);
fprintf('     Cohen d(paired)=%.3f\n\n', d_paired);

% --- Q3. 집단 x 조건 mixed ANOVA ---
fprintf('[Q3] 집단(ASD/TD) x 조건(FA/Hit) mixed ANOVA\n');
T = table(categorical(group), amp_err, amp_cor, ...
    'VariableNames', {'Group','FA','Hit'});
within = table(categorical({'FA';'Hit'}), 'VariableNames', {'Condition'});
rm = fitrm(T, 'FA-Hit ~ Group', 'WithinDesign', within);
ranovatbl = ranova(rm);                    % 조건 주효과 + 조건x집단 상호작용
btw = anova(rm);                           % 집단 주효과 (between)

disp('--- Between-subjects (집단 주효과) ---');
disp(btw);
disp('--- Within & Interaction (조건 주효과, 집단x조건) ---');
disp(ranovatbl);

% 상호작용을 ΔERN 집단비교로도 확인 (동등, 해석 쉬움)
d_asd = asd_amp_err - asd_amp_cor;
d_td  = td_amp_err  - td_amp_cor;
[~, p_int, ~, st_int] = ttest2(d_asd, d_td);
fprintf('\n[Q3-보조] ΔERN(FA-Hit) 집단비교 (independent t, 상호작용과 동등)\n');
fprintf('     ASD ΔERN = %.3f µV, TD ΔERN = %.3f µV\n', mean(d_asd), mean(d_td));
fprintf('     t(%d)=%.3f, p=%.4g\n', st_int.df, st_int.tstat, p_int);
fprintf('==========================\n');


%% ===================== Local functions =====================
function [err_mat, cor_mat, amp_err, amp_cor, times] = ...
         collect_group(file_list, proc_dir, roi_labels, bl_win, ern_win)
    err_mat=[]; cor_mat=[]; amp_err=[]; amp_cor=[]; times=[];
    for s = 1:numel(file_list)
        fn = char(string(file_list{s}));
        setfile = strrep(fn, '.mff', '_processed_data.set');
        try
            EEG = pop_loadset('filename', setfile, 'filepath', proc_dir);
        catch
            warning('로드 실패: %s', setfile); continue;
        end
        % ERP용 저역통과 필터 30 Hz (Wang 2020 traditional ERP: 0.1-30Hz).
        % 파이프라인 필터가 0.1-50Hz이므로 여기서 30Hz로 추가 저역통과.
        % (고주파 노이즈/근육 잔재 및 반응 실행 스파이크 완화). baseline correction 전에 적용.
        [~] = evalc('EEG = pop_eegfiltnew(EEG, [], 30);');
        EEG = pop_rmbase(EEG, bl_win);
        if isempty(times), times = EEG.times; end

        [idx_err, idx_cor] = get_condition_idx(EEG);
        if isempty(idx_err) || isempty(idx_cor)
            warning('조건 없음: %s', setfile); continue;
        end

        roi = zeros(1,numel(roi_labels));
        for c=1:numel(roi_labels)
            roi(c) = find(strcmpi({EEG.chanlocs.labels}, roi_labels{c}),1);
        end

        erp_err = mean(mean(EEG.data(roi,:,idx_err),3),1);  % [1 x time]
        erp_cor = mean(mean(EEG.data(roi,:,idx_cor),3),1);
        err_mat = [err_mat; erp_err];  %#ok<AGROW>
        cor_mat = [cor_mat; erp_cor];  %#ok<AGROW>

        % mean amplitude (ern_win 구간 평균)
        w = times >= ern_win(1) & times <= ern_win(2);
        amp_err = [amp_err; mean(erp_err(w))];  %#ok<AGROW>
        amp_cor = [amp_cor; mean(erp_cor(w))];  %#ok<AGROW>
    end
end

function [idx_err, idx_cor] = get_condition_idx(EEG)
    n_ep = EEG.trials;
    gonogo = nan(1,n_ep); accuracy = nan(1,n_ep);
    for e=1:n_ep
        lat = cell2mat(EEG.epoch(e).eventlatency);
        i0  = find(lat==0,1);
        if isempty(i0), continue; end
        gg = EEG.epoch(e).eventGoNogo;   if iscell(gg), gg=gg{i0}; end
        ac = EEG.epoch(e).eventAccuracy; if iscell(ac), ac=ac{i0}; end
        gonogo(e)=gg; accuracy(e)=ac;
    end
    idx_err = find(gonogo==2 & accuracy==0);
    idx_cor = find(gonogo==1 & accuracy==1);
end

function h = shaded(t, m, s, col)
    fill([t fliplr(t)], [m+s fliplr(m-s)], col, ...
        'FaceAlpha',0.2, 'EdgeColor','none', 'HandleVisibility','off');
    h = plot(t, m, 'Color', col, 'LineWidth', 1.8);
end

function decorate(ttl, ern_win)
    yl = get(gca,'ylim');
    % ERN 측정 구간 음영 표시
    patch([ern_win(1) ern_win(2) ern_win(2) ern_win(1)], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.9 0.6], 'FaceAlpha',0.25, 'EdgeColor','none', 'HandleVisibility','off');
    plot(get(gca,'xlim'), [0 0], 'k:', 'HandleVisibility','off');
    plot([0 0], yl, 'k:', 'HandleVisibility','off');
    xlabel('Time (ms)'); ylabel('\muV');
    title(ttl);
    % set(gca, 'YDir', 'reverse');   % negative를 위로
end