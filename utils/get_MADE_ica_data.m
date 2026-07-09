function [EEG, log, state] = get_MADE_ica_data(EEG, data_name, cfg, log, state, OVERLAP)
% GET_MADE_ICA_DATA  MADE STEP 8-11: ICA 준비 -> ICA 실행 -> ADJUST -> IC 제거
%
%   [EEG, log, state] = get_MADE_ica_data(EEG, data_name, cfg, log, state, OVERLAP)
%
%   get_MADE_filtered_data(STEP 1-7)를 통과한 EEG를 받아서 ICA 아티팩트
%   제거까지 수행한다. 패턴은 STEP 1-7 함수와 동일:
%     - checkpoint: 이미 ica_data가 있으면 load + log 복원 후 return
%     - all-bad-channels / all-bad-ICs 전멸 시 log.skip=true; return
%     - log는 데이터와 함께 저장 (.set -> EEG.etc.MADE_log, .mat -> log 변수)
%
%   입력 EEG는 STEP 7까지 끝난 상태(bad channel 제거 + reference 제거)라고 가정.

    output_dir = [cfg.output_location filesep 'ica_data(2)'];
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    [~, stem, ~] = fileparts(data_name);
    output_name = [output_dir filesep stem '_ica_data' cfg.output_format];

    % 이미 ICA까지 끝난 파일이 있으면 로드하고 반환 =====
    if ~OVERLAP && isfile(output_name)
        fprintf("\nPassing %s...(already ICA-cleaned)\n", data_name)
        if strcmp(cfg.output_format, '.set')
            EEG = pop_loadset('filename', [stem '_ica_data.set'], 'filepath', output_dir);
            if isfield(EEG, 'etc') && isfield(EEG.etc, 'MADE_log')
                log = EEG.etc.MADE_log; % log 복원
            else
                warning('MADE_log not found in %s. Using passed-in log.', [stem '_ica_data.set']);
            end
        elseif strcmp(cfg.output_format, '.mat')
            S = load([output_dir filesep stem '_ica_data.mat'], 'EEG', 'log');
            EEG = S.EEG;
            if isfield(S, 'log')
                log = S.log; % log 복원
            else
                warning('log not found in %s. Using passed-in log.', [stem '_ica_data.mat']);
            end
        end
        EEG = eeg_checkset(EEG);
        log.skip = false;
        return;
    end

    fprintf("\nRunning ICA for %s...\n", data_name)

    %% STEP 8: Prepare data for ICA
    EEG_copy = EEG;                 % 원본 복사
    EEG_copy = eeg_checkset(EEG_copy);

    % 복사본에 1Hz high-pass (ICA는 느린 drift에 취약하므로 강하게 필터)
    transband = 1;
    fl_cutoff = transband/2;
    fl_order  = 3.3 / (transband / EEG.srate);
    if mod(floor(fl_order),2) == 0
        fl_order = floor(fl_order);
    else
        fl_order = floor(fl_order)+1;
    end
    EEG_copy = pop_firws(EEG_copy, 'fcutoff', fl_cutoff, 'ftype', 'highpass', ...
        'wtype', 'hamming', 'forder', fl_order, 'minphase', 0);
    EEG_copy = eeg_checkset(EEG_copy);

    % [SPEED] ICA 학습용 사본만 다운샘플 (속도 ~2배).
    %   - 원본 EEG는 손대지 않으므로 최종 분석은 그대로 500Hz 유지.
    %   - EMG 탐지(20-40Hz)와 아티팩트 IC 분리에 필요한 대역은 250Hz(나이퀴스트 125Hz)로 충분.
    %   - 이미 1Hz highpass + (원본의 50Hz lowpass) 적용돼 있어 anti-aliasing 안전.
    ica_resample_rate = 250;   % Hz. 500->250. 필요 없으면 [] 로 두면 스킵.
    if ~isempty(ica_resample_rate) && EEG_copy.srate > ica_resample_rate
        EEG_copy = pop_resample(EEG_copy, ica_resample_rate);
        EEG_copy = eeg_checkset(EEG_copy);
    end

    % 1초 임시 epoch 생성 (dummy marker '999')
    EEG_copy = eeg_regepochs(EEG_copy, 'recurrence', 1, 'limits', [0 1], ...
        'rmbase', NaN, 'eventtype', '999');
    EEG_copy = eeg_checkset(EEG_copy);

    % 나쁜 epoch/채널 탐지 임계
    vol_thrs        = [-1000 1000];  % 전압 임계 (µV)
    emg_thrs        = [-100 30];     % 근전도(EMG) 파워 임계 (dB)
    emg_freqs_limit = [20 40];       % EMG 탐지 주파수 대역 (Hz)

    % 채널별로 아티팩트 epoch 비율 20% 초과하면 bad channel로 표시
    chanCounter = 1; ica_prep_badChans = [];
    numEpochs = EEG_copy.trials;

    for ch = 1:EEG_copy.nbchan
        % 전압 이상치 탐지
        EEG_copy = pop_eegthresh(EEG_copy, 1, ch, vol_thrs(1), vol_thrs(2), ...
            EEG_copy.xmin, EEG_copy.xmax, 0, 0);
        EEG_copy = eeg_checkset(EEG_copy);

        % 고주파(20-40Hz) EMG 탐지
        EEG_copy = pop_rejspec(EEG_copy, 1, 'elecrange', ch, 'method', 'fft', ...
            'threshold', emg_thrs, 'freqlimits', emg_freqs_limit, ...
            'eegplotplotallrej', 0, 'eegplotreject', 0);
        EEG_copy = eeg_checkset(EEG_copy);

        EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
        artifacted_epochs = EEG_copy.reject.rejglobal;

        if sum(artifacted_epochs) > (numEpochs*20/100)
            ica_prep_badChans(chanCounter) = ch; %#ok<AGROW>
            chanCounter = chanCounter + 1;
        end
    end

    % bad channel 라벨 기록 (반드시 pop_select로 제거하기 '전'에)
    % 이 시점 EEG_copy는 아직 채널 제거 전이라 인덱스가 정확히 맞음
    % 인덱스는 채널 구성에 종속돼 헷갈리므로 라벨로만 저장
    if numel(ica_prep_badChans) == 0
        log.ica_preparation_bad_channels = 'none';
    else
        log.ica_preparation_bad_channels = strjoin({EEG_copy.chanlocs(ica_prep_badChans).labels}, ',');
    end

    % 전멸 판정: (거의) 모든 채널이 bad
    all_bad_channels = 0;
    if numel(ica_prep_badChans) == EEG.nbchan || numel(ica_prep_badChans)+1 == EEG.nbchan
        all_bad_channels = 1;
        warning(['No usable data (all bad channels in ICA prep) for datafile ', data_name]);
    else
        % bad channel 제거 (복사본)
        EEG_copy = pop_select(EEG_copy, 'nochannel', ica_prep_badChans);
        EEG_copy = eeg_checkset(EEG_copy);
    end

    % 전멸이면: log 채우고 데이터 저장 후 skip
    if all_bad_channels == 1
        log.length_ica_data                        = 0;
        log.total_ICs                              = 0;
        log.ICs_removed                            = '0';
        log.total_epochs_before_artifact_rejection = 0;
        log.total_epochs_after_artifact_rejection  = 0;
        log.total_channels_interpolated            = 0;
        log.ica_data_rank                          = 0;
        log.adjust_performed                       = 'no_ica';   % ICA 자체를 못 함
        log.skip = true;
        save_all_bad(EEG, log, cfg, output_dir, stem, '_no_usable_data_all_bad_channels');
        return;
    end

    % 전체 채널 기준으로 아티팩트 epoch 찾아서 ICA 전에 제거
    EEG_copy = pop_eegthresh(EEG_copy, 1, 1:EEG_copy.nbchan, vol_thrs(1), vol_thrs(2), ...
        EEG_copy.xmin, EEG_copy.xmax, 0, 0);
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = pop_rejspec(EEG_copy, 1, 'elecrange', 1:EEG_copy.nbchan, 'method', 'fft', ...
        'threshold', emg_thrs, 'freqlimits', emg_freqs_limit, ...
        'eegplotplotallrej', 0, 'eegplotreject', 0);
    EEG_copy = eeg_checkset(EEG_copy);
    EEG_copy = eeg_rejsuperpose(EEG_copy, 1, 1, 1, 1, 1, 1, 1, 1);
    reject_artifacted_epochs = EEG_copy.reject.rejglobal;
    EEG_copy = pop_rejepoch(EEG_copy, reject_artifacted_epochs, 0);

    %% STEP 9: Run ICA
    log.length_ica_data = EEG_copy.trials;   % ICA에 들어간 epoch 수
    EEG_copy = eeg_checkset(EEG_copy);

    % [RANK] 데이터의 실제 rank를 계산해서 PCA 차원으로 명시
    %   reference 제거(예: E65)로 rank가 채널수보다 작으면, PCA 미지정 시 runica가
    %   IC를 채널수보다 적게 만들어 icaweights가 직사각이 되고 -> STEP 10 ADJUST가 스킵됨
    %   'pca'로 IC 수를 rank에 맞추면 icaweights가 rank×rank 정사각이 되어 ADJUST가 항상 작동
    %   평균 제거 후 계산해 EEG 채널 간 높은 상관에 의한 rank 과대추정을 줄임
    tmp_for_rank = EEG_copy.data(:, :);
    tmp_for_rank = double(tmp_for_rank) - mean(double(tmp_for_rank), 2);
    data_rank = rank(tmp_for_rank);
    log.ica_data_rank = data_rank;   % 참고용 기록
    clear tmp_for_rank;

    EEG_copy = pop_runica(EEG_copy, 'icatype', 'runica', 'extended', 1, ...
        'pca', data_rank, 'stop', 1E-7, 'interupt', 'off');

    % 복사본에서 구한 ICA 가중치를 원본에 이식
    ICA_WINV     = EEG_copy.icawinv;
    ICA_SPHERE   = EEG_copy.icasphere;
    ICA_WEIGHTS  = EEG_copy.icaweights;
    ICA_CHANSIND = EEG_copy.icachansind;

    % 복사본에서 제거한 bad channel을 원본에서도 동일하게 제거 (채널 구성 일치)
    EEG = eeg_checkset(EEG);
    EEG = pop_select(EEG, 'nochannel', ica_prep_badChans);

    EEG.icawinv     = ICA_WINV;
    EEG.icasphere   = ICA_SPHERE;
    EEG.icaweights  = ICA_WEIGHTS;
    EEG.icachansind = ICA_CHANSIND;
    EEG = eeg_checkset(EEG);

    %% STEP 10: ADJUST로 나쁜 IC 자동 탐지
    badICs = [];
    EEG_copy = EEG;
    EEG_copy = eeg_regepochs(EEG_copy, 'recurrence', 1, 'limits', [0 1], ...
        'rmbase', NaN, 'eventtype', '999');
    EEG_copy = eeg_checkset(EEG_copy);

    % [RANK-DEFICIENT 대응] 원본 MADE는 icaweights가 정사각(IC수=채널수)일 때만
    % ADJUST를 실행. 그러나 reference 제거로 rank가 채널수보다 작으면 pop_runica가
    % rank를 올바르게 감지해 IC를 rank만큼만 생성하므로(=ghost IC 없음, 건강한 상태)
    % icaweights는 항상 IC수×채널수 직사각이 됨. 이 경우에도 ADJUST는 icawinv(채널×IC)
    % 기반으로 IC 수를 세어 정상 동작함. 따라서 정사각 대신 "IC가 2개 이상"인지로 판단함.
    if size(EEG_copy.icawinv, 2) >= 2
        adjust_report = [output_dir filesep stem '_adjust_report'];

        % [ICAACT 직접 계산] ADJUST의 icaact 재계산 코드(adjusted_ADJUST 154번)는
        % IC수=채널수를 가정해 size(EEG.data,1)로 reshape하므로, pca로 rank 축소된
        % 경우(IC수<채널수) 요소 개수 불일치로 reshape 에러가 남.
        % 따라서 icachansind(ICA에 실제 쓰인 채널)만 골라 icaact를 올바른 IC 차원으로
        % 미리 계산해 채워두면, ADJUST가 재계산을 건너뛰어 에러를 피함
        icaact_2d = (EEG_copy.icaweights * EEG_copy.icasphere) ...
                    * EEG_copy.data(EEG_copy.icachansind, :);
        EEG_copy.icaact = reshape(icaact_2d, size(icaact_2d,1), ...
                                  EEG_copy.pnts, EEG_copy.trials);

        % [TRY-CATCH] ADJUST 내부가 실패해도(EM 'Negative Discriminant' 등 데이터 분포
        % 문제, 또는 기타 내부 에러) 파이프라인이 멈추지 않도록 감쌈. 실패한 subject는
        % IC를 제거하지 않고(badICs=[]) adjust_performed에 사유를 기록해 나중에 수동 검수.
        try
            badICs = adjusted_ADJUST(EEG_copy, adjust_report);
            log.adjust_performed = 'yes';   % ADJUST 정상 실행
        catch ME
            badICs = [];
            reason = strtrim(ME.message);
            reason = reason(1:min(60, numel(reason)));   % 메시지 앞부분만
            log.adjust_performed = ['failed: ' reason];
            warning('ADJUST failed for %s (%s). No ICs removed; flagged for manual review.', ...
                data_name, reason);
        end
        close all;
    else
        % IC가 2개 미만이면 ADJUST 불가 (극히 드묾)
        log.adjust_performed = 'skipped_no_ica';
        warning(['ICA produced fewer than 2 components. ADJUST skipped for ', data_name]);
    end

    % 나쁜 IC 표시
    for ic = 1:length(badICs)
        EEG.reject.gcompreject(1, badICs(ic)) = 1;
        EEG = eeg_checkset(EEG);
    end
    log.total_ICs = size(EEG.icasphere, 1);
    if numel(badICs) == 0
        log.ICs_removed = '0';
    else
        log.ICs_removed = num2str(double(badICs));
    end

    % interim 저장 (ICA 가중치 이식된 상태, IC 제거 전) — checkpoint용은 아래 최종 저장
    % (원본은 여기서 _ica_data를 저장하지만, 본 함수는 IC 제거 후 최종본을 저장하여
    %  checkpoint로 삼는다. IC 제거 전 상태가 별도로 필요하면 여기서 저장 추가 가능.)
    % IC 제거 후 데이터 이상하게 생겼는지 확인하기 위해 랜덤으로 몇개만 하기

    %% STEP 11: 나쁜 IC 실제 제거
    all_bad_ICs = 0;
    ICs2remove = find(EEG.reject.gcompreject);

    if numel(ICs2remove) == log.total_ICs
        % 모든 IC가 나쁨 -> 쓸 데이터 없음
        all_bad_ICs = 1;
        warning(['No usable data (all bad ICs) for datafile ', data_name]);
        log.total_epochs_before_artifact_rejection = 0;
        log.total_epochs_after_artifact_rejection  = 0;
        log.total_channels_interpolated            = 0;
        log.skip = true;
        save_all_bad(EEG, log, cfg, output_dir, stem, '_no_usable_data_all_bad_ICs');
        return;
    else
        EEG = eeg_checkset(EEG);
        EEG = pop_subcomp(EEG, ICs2remove, 0);   % 나쁜 IC 제거
    end

    %% 최종 저장 (checkpoint) — IC 제거 완료된 상태
    log.skip = false;
    if cfg.save_interim_result == 1
        if strcmp(cfg.output_format, '.set')
            EEG = eeg_checkset(EEG);
            EEG.etc.MADE_log = log;
            EEG = pop_editset(EEG, 'setname', [stem '_ica_data']);
            EEG = pop_saveset(EEG, 'filename', [stem '_ica_data.set'], 'filepath', output_dir);
        elseif strcmp(cfg.output_format, '.mat')
            save([output_dir filesep stem '_ica_data.mat'], 'EEG', 'log');
        end
    end
end


function save_all_bad(EEG, log, cfg, output_dir, stem, suffix)
% 전멸(all bad) 상태의 데이터 + log를 저장하는 헬퍼
    if strcmp(cfg.output_format, '.set')
        EEG = eeg_checkset(EEG);
        EEG.etc.MADE_log = log;
        EEG = pop_editset(EEG, 'setname', [stem suffix]);
        EEG = pop_saveset(EEG, 'filename', [stem suffix '.set'], 'filepath', output_dir);
    elseif strcmp(cfg.output_format, '.mat')
        save([output_dir filesep stem suffix '.mat'], 'EEG', 'log');
    end
end
