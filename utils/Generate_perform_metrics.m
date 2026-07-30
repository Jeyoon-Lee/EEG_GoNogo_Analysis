function generate_perform_metrics(filtered_dir, report_path)
% GENERATE_PERFORM_METRICS  Go/No-Go 수행지표를 계산해 CSV로 저장
%
%   generate_perform_metrics(filtered_dir, report_path)
%
%   Add_block_event로 block/trial 태깅이 끝난 filtered_data(.set)를 읽어서
%   Go/No-Go 수행지표를 계산하고 CSV에 한 줄씩 append
%
%   입력:
%     filtered_dir : filtered_data(.set)가 있는 폴더
%                    예: 'ERP/Data/T2/filtered_data'
%     report_path  : 결과 CSV 경로
%                    예: 'ERP/Data/T2/performance_metrics.csv'
%
%   trial 판정 (trial_global 로 이벤트를 묶어서):
%     Go(cel=101)   + resp 있음 -> Hit,  RT = resp_lat - stm_lat
%     Go(cel=101)   + resp 없음 -> OE  (omission error)
%     NoGo(cel=102) + resp 있음 -> FA  (false alarm), FA_RT 계산
%     NoGo(cel=102) + resp 없음 -> CR  (correct rejection)
%
%   RT valid 범위: Go Hit 한정 100~1500ms. (FA_RT는 raw)

    %% ===== Parameters =====
    cel_marker  = 'stm+';      % cel(101/102)을 읽을 마커 (한 trial 내 모두 동일)
    stim_marker = 'stm+';      % RT 기준 자극 마커
    resp_marker = 'resp';
    go_code     = '101';
    nogo_code   = '102';
    rt_valid    = [100 1500];  % Go Hit RT 유효범위 (ms)

    files = dir(fullfile(char(filtered_dir), '*_filtered_data.set'));
    fprintf('%d개 파일 처리 시작\n', numel(files));

    skipped = {};

    for i = 1:numel(files)
        fname = files(i).name;
        fprintf('\n(%d/%d) %s\n', i, numel(files), fname);

        EEG = pop_loadset('filename', fname, 'filepath', char(filtered_dir));

        % block/trial 태깅 여부 확인
        if ~isfield(EEG.event, 'trial_global')
            skipped{end+1} = sprintf('%s: trial_global 없음 (태깅 안 됨)', fname); %#ok<AGROW>
            warning('  trial_global 필드 없음. 건너뜀.');
            continue;
        end

        types_all = {EEG.event.type};
        tg        = [EEG.event.trial_global];
        utrials   = unique(tg(tg > 0));

        % ---- 지표 누적 ----
        Hit=0; OE=0; FA=0; CR=0;
        RT_hit_all = [];   % 모든 Hit RT (raw)
        RT_fa_all  = [];   % 모든 FA RT (raw)

        for t = utrials(:)'
            idx = find(tg == t);
            tt  = types_all(idx);

            % 이 trial의 조건 (cel_marker에서 읽기)
            ci = find(strcmp(tt, cel_marker), 1);
            if isempty(ci), continue; end   % 자극 마커 없는 이상 trial 건너뜀
            cel = EEG.event(idx(ci)).mffkey_cel;
            if isnumeric(cel), cel = num2str(cel); end

            % resp 유무 및 RT
            ri = find(strcmp(tt, resp_marker), 1);
            si = find(strcmp(tt, stim_marker), 1);
            has_resp = ~isempty(ri);

            rt = NaN;
            if has_resp && ~isempty(si)
                rt = (EEG.event(idx(ri)).latency - EEG.event(idx(si)).latency) / EEG.srate * 1000;
            end

            if strcmp(cel, go_code)             % Go
                if has_resp
                    Hit = Hit + 1;
                    RT_hit_all(end+1) = rt; %#ok<AGROW>
                else
                    OE = OE + 1;
                end
            elseif strcmp(cel, nogo_code)       % NoGo
                if has_resp
                    FA = FA + 1;
                    RT_fa_all(end+1) = rt; %#ok<AGROW>
                else
                    CR = CR + 1;
                end
            end
        end

        n_Go   = Hit + OE;
        n_NoGo = FA + CR;
        total_trials = numel(utrials);
        n_Blocks = max([EEG.event.block]);

        % ---- RT 통계 ----
        % raw (모든 Hit)
        RT_mean = mean_or_nan(RT_hit_all);
        RT_sd   = std_or_nan(RT_hit_all);
        % valid (100~1500ms Hit만)
        valid_mask = RT_hit_all >= rt_valid(1) & RT_hit_all <= rt_valid(2);
        RT_hit_valid = RT_hit_all(valid_mask);
        RT_mean_valid = mean_or_nan(RT_hit_valid);
        RT_sd_valid   = std_or_nan(RT_hit_valid);
        n_RT_invalid  = sum(~valid_mask);   % 범위 밖 Hit 개수
        % FA RT (raw)
        FA_RT = mean_or_nan(RT_fa_all);

        % 비율
        Hit_rate = safe_div(Hit, n_Go);
        FA_rate  = safe_div(FA, n_NoGo);

        % ID / phase 파싱
        id_tok    = regexp(fname, '[A-Z]{2}\d{3}', 'match', 'once');
        phase_tok = regexp(fname, '_(T\d)_', 'tokens', 'once');
        ID    = string(id_tok);
        phase = string(phase_tok);

        % CSV 한 줄 append
        row = table( ...
            ID, phase, n_Blocks, total_trials, n_Go, n_NoGo, ...
            Hit, OE, FA, CR, Hit_rate, FA_rate, ...
            RT_mean, RT_sd, RT_mean_valid, RT_sd_valid, n_RT_invalid, FA_RT, ...
            'VariableNames', {'ID','phase','n_Blocks','total_trials','n_Go','n_NoGo', ...
                'Hit','OE','FA','CR','Hit_rate','FA_rate', ...
                'RT_mean','RT_sd','RT_mean_valid','RT_sd_valid','n_RT_invalid','FA_RT'});

        if isfile(char(report_path))
            writetable(row, char(report_path), 'WriteMode', 'append');
        else
            writetable(row, char(report_path));
        end

        fprintf('  Hit=%d OE=%d FA=%d CR=%d | RT_valid=%.0f±%.0f (invalid %d) | FA_RT=%.0f\n', ...
            Hit, OE, FA, CR, RT_mean_valid, RT_sd_valid, n_RT_invalid, FA_RT);
    end

    %% 확인 필요 subject 출력
    fprintf('\n===== 확인 필요 subject =====\n');
    if isempty(skipped)
        fprintf('없음\n');
    else
        for s = 1:numel(skipped), fprintf('%s\n', skipped{s}); end
    end
    fprintf('\nDone\n');
end

%% ===== Helpers =====
function m = mean_or_nan(x)
    if isempty(x), m = NaN; else, m = mean(x); end
end

function s = std_or_nan(x)
    if numel(x) < 2, s = NaN; else, s = std(x); end
end

function r = safe_div(a, b)
    if b == 0, r = NaN; else, r = a / b; end
end