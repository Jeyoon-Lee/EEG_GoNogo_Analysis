function [EEG, log, state] = get_MADE_processed_data(EEG, data_name, cfg, log, state, OVERLAP)
% GET_MADE_PROCESSED_DATA  MADE STEP 12-16: epoching -> baseline -> artifact
%   rejection(+epoch interpolation) -> channel interpolation -> rereference
%
%   [EEG, log, state] = get_MADE_processed_data(EEG, data_name, cfg, log, state, OVERLAP)
%
%   ICA(STEP 8-11)까지 끝난 EEG를 받아서 분석용 epoch 데이터를 만든다.
%
%   설계 메모:
%   - epoch은 response('resp') 기준으로 자른다. 각 epoch에 GoNogo/Accuracy/RT 라벨을
%     심어두므로, 나중에 NoGo-FA(ERN) / Go-Hit(CRN) 등 조건은 분석 단계에서 뽑는다.
%     (미리 조건별로 파일을 나누지 않음 → 한 파일에서 모든 조건 추출 가능)
%   - frontal_channels 기반 epoch 사전 제거는 사용하지 않는다. ADJUST가 눈 아티팩트
%     IC를 이미 제거하므로 중복이며 과도한 epoch 손실을 유발할 수 있다.
%   - 파라미터(epoch 창, threshold 등)를 바꿔가며 여러 번 돌리는 단계이므로,
%     결과는 processed_tag 별 폴더에 저장해 조건별로 구분한다.
%    stim_marker는 stim+, resp_marker는 resp 로 고정!

    % ===== 저장 경로 (조건별 태그로 구분) =====
    processed_tag = string(cfg.task_event_markers{1}) + "_" + strjoin(string(cfg.task_epoch_length * 1000), "_"); % 이 과정 때문에 epoch를 자르는 marker 기준은 하나여야 함
    out_sub = ['processed_data_' char(processed_tag)];
    output_dir = [cfg.output_location filesep out_sub];
    if ~isfolder(output_dir), mkdir(output_dir); end

    [~, stem, ~] = fileparts(data_name);
    output_name = [output_dir filesep stem '_processed_data' cfg.output_format];

    % ===== Checkpoint =====
    if ~OVERLAP && isfile(output_name)
        fprintf("\nPassing %s...(already processed)\n", data_name)
        if strcmp(cfg.output_format, '.set')
            EEG = pop_loadset('filename', [stem '_processed_data.set'], 'filepath', output_dir);
            if isfield(EEG, 'etc') && isfield(EEG.etc, 'MADE_log')
                log = EEG.etc.MADE_log;
            end
        elseif strcmp(cfg.output_format, '.mat')
            S = load([output_dir filesep stem '_processed_data.mat'], 'EEG', 'log');
            EEG = S.EEG;
            if isfield(S, 'log'), log = S.log; end
        end
        EEG = eeg_checkset(EEG);
        log.skip = false;
        return;
    end

    fprintf("\nProcessing (STEP 12-16) %s...\n", data_name)

    %% STEP 12: 이벤트 라벨링 + epoching
    % ---- 12a. 이벤트에 GoNogo/Accuracy/RT 라벨 심기 (trial_global 기반) ----
    [EEG.event.StimType]  = deal(NaN);
    [EEG.event.GoNogo]    = deal(NaN);
    [EEG.event.Responded] = deal(NaN);
    [EEG.event.Accuracy]  = deal(NaN);
    [EEG.event.RT]        = deal(NaN);

    types_all = {EEG.event.type};
    tg        = [EEG.event.trial_global];
    utrials   = unique(tg(tg > 0));

    for t = utrials(:)'
        idx = find(tg == t);
        tt  = types_all(idx);

        si = find(strcmp(tt, 'stm+'), 1);   % 자극 마커
        ri = find(strcmp(tt, 'resp'), 1);   % 응답 마커 (첫 번째)
        if isempty(si), continue; end
        si_ev = idx(si);

        cel = EEG.event(si_ev).mffkey_cel;
        if isnumeric(cel), cel = num2str(cel); end
        is_go   = strcmp(cel, '101');
        is_nogo = strcmp(cel, '102');

        % stim 마커 라벨
        EEG.event(si_ev).StimType = 1;
        if is_go,   EEG.event(si_ev).GoNogo = 1; end
        if is_nogo, EEG.event(si_ev).GoNogo = 2; end

        if ~isempty(ri)
            ri_ev = idx(ri);
            % resp 마커 라벨
            EEG.event(ri_ev).StimType = 2;
            if is_go,   EEG.event(ri_ev).GoNogo = 1; end
            if is_nogo, EEG.event(ri_ev).GoNogo = 2; end

            EEG.event(si_ev).Responded = 1;
            EEG.event(ri_ev).Responded = 1;

            % Accuracy: Go+resp=Hit(1), NoGo+resp=FA=error(0)
            if is_go
                acc = 1;
            elseif is_nogo
                acc = 0;
            else
                acc = NaN;
            end
            EEG.event(si_ev).Accuracy = acc;
            EEG.event(ri_ev).Accuracy = acc;

            % RT (ms)
            rt = (EEG.event(ri_ev).latency - EEG.event(si_ev).latency) * (1000/EEG.srate);
            EEG.event(si_ev).RT = rt;
            EEG.event(ri_ev).RT = rt;
        else
            % 응답 없음
            EEG.event(si_ev).Responded = 0;
            % Accuracy: Go+noresp=miss(0), NoGo+noresp=CR(1)
            if is_go
                EEG.event(si_ev).Accuracy = 0;
            elseif is_nogo
                EEG.event(si_ev).Accuracy = 1;
            end
        end
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');

    % ---- 12b. response 기준 epoch 자르기 ----
    EEG = eeg_checkset(EEG);
    EEG = pop_epoch(EEG, cfg.task_event_markers, cfg.task_epoch_length, 'epochinfo', 'yes');
    log.total_epochs_before_artifact_rejection = EEG.trials;

    %% STEP 13: baseline 제거
    if cfg.remove_baseline == 1
        EEG = eeg_checkset(EEG);
        EEG = pop_rmbase(EEG, cfg.baseline_window);
    end

    %% STEP 14: Artifact rejection + epoch-level interpolation
    % (frontal_channels 사전 제거 없음: ADJUST가 눈 아티팩트 처리)
    all_bad_epochs = 0;
    if cfg.voltthres_rejection == 1
        if cfg.interp_epoch == 1
            % 각 채널에서 전압 초과 epoch 탐지 (제거하지 않고 표시만)
            % [속도] 채널별 eeg_checkset 제거(불필요), pop_eegthresh 로그는 evalc로 억제.
            badChans = zeros(EEG.nbchan, EEG.trials);
            for ch = 1:EEG.nbchan
                [~] = evalc(['EEG = pop_eegthresh(EEG, 1, ch, cfg.volt_threshold(1), ' ...
                    'cfg.volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);']);
                EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                badChans(ch, :) = EEG.reject.rejglobal;
            end

            % 먼저 "버릴 epoch"을 판정한다: 한 epoch에서 나쁜 채널이 전체의 10%를
            % 넘으면 그 epoch은 어차피 제거 대상이다. 이런 epoch은 보간을 시도하지 않는다.
            % (거의 모든 채널이 나쁜 epoch을 보간하려 하면 좋은 채널이 없어 구면스플라인
            %  계산이 NaN이 되어 legendre에서 에러가 나기 때문)
            reject_thresh = round((10/100) * EEG.nbchan);
            badepoch = zeros(1, EEG.trials);
            for ei = 1:EEG.trials
                if sum(badChans(:, ei)) > reject_thresh
                    badepoch(ei) = 1;
                end
            end
            badepoch = logical(badepoch);

            % epoch별로 나쁜 채널만 보간 (단, 버릴 epoch은 스킵)
            % pop_selectevent로 epoch을 정식으로 떼어내 보간(안정적). 로그는 evalc로 억제.
            tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
            for e = 1:EEG.trials
                badChanNum = find(badChans(:, e) == 1);
                if isempty(badChanNum) || badepoch(e)
                    % 보간할 것이 없거나, 어차피 버릴 epoch이면 원본 그대로 두고 넘어감
                    tmpData(:, :, e) = EEG.data(:, :, e);
                    continue;
                end
                [~, EEGe] = evalc(['pop_selectevent(EEG, ''epoch'', e, ' ...
                    '''deleteevents'', ''off'', ''deleteepochs'', ''on'', ''invertepochs'', ''off'')']);
                [~, EEGe_interp] = evalc('eeg_interp(EEGe, badChanNum)');
                tmpData(:, :, e) = EEGe_interp.data;
            end
            EEG.data = tmpData;

            % badepoch은 위(보간 전)에서 이미 계산됨. 그 판정으로 epoch 제거.
            if sum(badepoch) == EEG.trials || sum(badepoch)+1 == EEG.trials
                all_bad_epochs = 1;
            else
                EEG = pop_rejepoch(EEG, badepoch, 0);
                EEG = eeg_checkset(EEG);
            end
        else
            % epoch interpolation 없이 전압 초과 epoch 그냥 제거
            EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, cfg.volt_threshold(1), cfg.volt_threshold(2), ...
                EEG.xmin, EEG.xmax, 0, 0);
            EEG = eeg_checkset(EEG);
            EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
            if sum(EEG.reject.rejthresh) == EEG.trials || sum(EEG.reject.rejthresh)+1 == EEG.trials
                all_bad_epochs = 1;
            else
                EEG = pop_rejepoch(EEG, EEG.reject.rejthresh, 0);
                EEG = eeg_checkset(EEG);
            end
        end
    end

    % 전멸 처리
    if all_bad_epochs == 1
        warning(['No usable data (all bad epochs) for ', data_name]);
        log.total_epochs_after_artifact_rejection = 0;
        log.total_channels_interpolated           = 0;
        log.n_epoch_go_hit  = 0;
        log.n_epoch_go_miss = 0;
        log.n_epoch_nogo_fa = 0;
        log.n_epoch_nogo_cr = 0;
        log.n_epoch_uncat   = 0;
        log.skip = true;
        save_processed(EEG, log, cfg, output_dir, stem, '_no_usable_data_all_bad_epochs');
        return;
    end
    log.total_epochs_after_artifact_rejection = EEG.trials;

    % ---- 조건별 살아남은 epoch 수 집계 (ERN/CRN 스크리닝용) ----
    % 각 epoch에서 기준 이벤트(latency==0, 즉 자른 그 resp)의 GoNogo/Accuracy 라벨로 분류.
    %   GoNogo: 1=Go, 2=NoGo / Accuracy: Go 1=Hit,0=miss / NoGo 0=FA(error),1=CR
    n_go_hit = 0; n_go_miss = 0; n_nogo_fa = 0; n_nogo_cr = 0; n_uncat = 0;
    for e = 1:EEG.trials
        lat = cell2mat(EEG.epoch(e).eventlatency);
        i0  = find(lat == 0, 1);                    % 기준 이벤트
        if isempty(i0), n_uncat = n_uncat + 1; continue; end
        gg = EEG.epoch(e).eventGoNogo;   if iscell(gg), gg = gg{i0}; end
        ac = EEG.epoch(e).eventAccuracy; if iscell(ac), ac = ac{i0}; end
        if     gg==1 && ac==1, n_go_hit  = n_go_hit  + 1;   % Go-Hit (CRN)
        elseif gg==1 && ac==0, n_go_miss = n_go_miss + 1;   % Go-miss
        elseif gg==2 && ac==0, n_nogo_fa = n_nogo_fa + 1;   % NoGo-FA (ERN)
        elseif gg==2 && ac==1, n_nogo_cr = n_nogo_cr + 1;   % NoGo-CR
        else,                  n_uncat   = n_uncat   + 1;
        end
    end
    log.n_epoch_go_hit  = n_go_hit;    % Go-Hit  (CRN 대상)
    log.n_epoch_go_miss = n_go_miss;   % Go-miss
    log.n_epoch_nogo_fa = n_nogo_fa;   % NoGo-FA (ERN 대상) ★핵심
    log.n_epoch_nogo_cr = n_nogo_cr;   % NoGo-CR
    log.n_epoch_uncat   = n_uncat;     % 미분류 (라벨 없는 기준이벤트 등)

    %% STEP 15: 제거된 bad channel 다시 보간
    % state.channels_analysed에 의존하지 않고, 원본 채널 위치(.ced)에서 목표 채널을
    % 직접 구성해 보간한다. checkpoint로 load된 경우에도 안전.
    % .ced에는 실제 전극(E1~E64, type=EEG) 외에 reference(E65, type=EEG)와
    % fiducial(E66~E68, type=FID)이 섞여 있다. 실제 두피 전극만 목표로 해야 하므로
    % type=='EEG' 이면서 reference가 아닌 채널(=E1~E64, 64개)만 남긴다.
    if cfg.interp_channels == 1
        n_before_interp = EEG.nbchan;   % 보간 전 채널 수 (몇 개 복원되는지 계산용)

        full_chanlocs = readlocs(char(cfg.channel_locations));   % 원본 전체 채널 위치

        ref_label = 'E65';
        if isfield(cfg, 'reference_label') && ~isempty(cfg.reference_label)
            ref_label = cfg.reference_label;
        end

        % type이 'EEG'인 것만 (fiducial 'FID' 등 제외). type 필드가 없을 경우 대비.
        if isfield(full_chanlocs, 'type')
            is_eeg = strcmpi({full_chanlocs.type}, 'EEG');
        else
            is_eeg = true(1, numel(full_chanlocs));
        end
        is_ref = strcmpi({full_chanlocs.labels}, ref_label);   % reference 제외
        keep = is_eeg & ~is_ref;
        target_chanlocs = full_chanlocs(keep);   % E1~E64 (64채널) 목표 위치

        EEG = eeg_interp(EEG, target_chanlocs);   % 현재 없는 채널을 보간해 복원
        EEG = eeg_checkset(EEG);

        % 실제로 보간(복원)된 채널 수 = 보간 후 - 보간 전
        log.total_channels_interpolated = EEG.nbchan - n_before_interp;
    else
        log.total_channels_interpolated = 0;
    end

    %% STEP 16: Rereference
    if cfg.rerefer_data == 1
        EEG = eeg_checkset(EEG);
        if iscell(cfg.reref)
            reref_idx = zeros(1, length(cfg.reref));
            for rr = 1:length(cfg.reref)
                reref_idx(rr) = find(strcmp({EEG.chanlocs.labels}, cfg.reref{rr}));
            end
            EEG = pop_reref(EEG, reref_idx);
        else
            EEG = pop_reref(EEG, cfg.reref);   % [] 이면 average reference
        end
    end

    %% 최종 저장 (checkpoint)
    log.skip = false;
    if strcmp(cfg.output_format, '.set')
        EEG = eeg_checkset(EEG);
        EEG.etc.MADE_log = log;
        EEG = pop_editset(EEG, 'setname', [stem '_processed_data']);
        EEG = pop_saveset(EEG, 'filename', [stem '_processed_data.set'], 'filepath', output_dir);
    elseif strcmp(cfg.output_format, '.mat')
        save([output_dir filesep stem '_processed_data.mat'], 'EEG', 'log');
    end
end


function save_processed(EEG, log, cfg, output_dir, stem, suffix)
    if strcmp(cfg.output_format, '.set')
        EEG = eeg_checkset(EEG);
        EEG.etc.MADE_log = log;
        EEG = pop_editset(EEG, 'setname', [stem suffix]);
        EEG = pop_saveset(EEG, 'filename', [stem suffix '.set'], 'filepath', output_dir);
    elseif strcmp(cfg.output_format, '.mat')
        save([output_dir filesep stem suffix '.mat'], 'EEG', 'log');
    end
end
