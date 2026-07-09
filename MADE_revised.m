clear;
clc;
eeglab;
addpath(genpath(pwd));
overwrite = false;


% 추후 state 제거하기
% reference channel E65로 고정
%% %%%%%%%%%%%%%%%%%%%%%%%%    PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = struct();
cfg.rawdata_location        = 'ERP/Raw_Data/T2';
cfg.output_location         = 'ERP/Data/T2';
cfg.channel_locations       = 'ERP/Raw_Data_Info/download65channels.ced';
cfg.report_name             = 'MADE_report_260704.csv';

cfg.adjust_time_offset      = true; % need to be checked
    cfg.filter_timeoffset   = 0; % ms
    cfg.stimulus_timeoffset = 41.78; % ms
    cfg.response_timeoffset = 0; % ms
    cfg.stimulus_markers    = {'fix+', 'stm+', 'bgin'};
    cfg.response_markers    = {'resp'};

cfg.down_sample             = false;
cfg.sampling_rate           = 500; % Hz
cfg.delete_outerlayer       = false;
    cfg.outerlayer_channel  = {};

cfg.highpass                = 0.3; % Hz
cfg.lowpass                 = 50; % Hz

cfg.epoch_data              = true;
cfg.task_event_markers      = {'stm+'}; % {'resp'} or {'stim+'} % 0 sec 기준점
cfg.task_epoch_length       = [-0.2 0.8]; % [-0.4 0.6] for resp, [-0.2 0.8] for stim+ % second

cfg.rest_epoch_length       = 0;
cfg.overlap_epoch           = false;
cfg.dummy_events            = {'999'};

cfg.remove_baseline         = false;
cfg.baseline_window         = []; % changed from [-400 -200]

cfg.voltthres_rejection     = true;
cfg.volt_threshold          = [-125 125];

cfg.interp_epoch            = true;
cfg.interp_channels         = true;

cfg.rerefer_data            = true;
cfg.reref                   = []; % [] means average rereference

cfg.save_interim_result     = true;
cfg.output_format           = '.set';  % '.set' or '.mat' (Matlab data structure)

% Preprocessing

datafile_names = dir([cfg.rawdata_location filesep '*.mff']);

% 실행 시작 시 기존 리포트 삭제 (매번 깨끗하게 새로 append)
report_path = [cfg.output_location filesep cfg.report_name];
if isfile(report_path)
    delete(report_path);
end

for subj = 1:length(datafile_names)
    cur_dname = datafile_names(subj).name;

    log = init_log();
    clc;
    fprintf("\n\n\n(%d/%d)", subj, length(datafile_names))
    
    %% get_MADE_filtered_data (-STEP 7)
    % import -> channel location -> time offset adjustment -> resampling
    % -> delete outerlayer channels -> filtering -> FASTER bad channel detection
    % output save at output_location/filtered_data
    [EEG, log, state] = get_MADE_filtered_data(cur_dname, cfg, log, overwrite);

    if log.skip
        append_report_row(log, cur_dname, cfg);
        continue;
    end
    
    EEG = Add_block_event(EEG, cur_dname, cfg);

    %% get_MADE_ica_data (STEP 8-11)
    % ICA 준비 -> ICA 실행 _> ADJUST -> IC 제거
    % output save at ouput_location/ica_data
    % 반드시 함수 내부 output_dir 확인할 것
    [EEG,log, state] = get_MADE_ica_data(EEG, cur_dname, cfg, log, state, overwrite);

    if log.skip
        append_report_row(log, cur_dname, cfg);
        continue;
    end

    % get_MADE_processed_data (STEP 12-)
    % 기준 marker stm+, resp로 고정
    [EEG, log, state] = get_MADE_processed_data(EEG, cur_dname, cfg, log, state, overwrite);
    if log.skip
        append_report_row(log, cur_dname, cfg);
        continue;
    end
    % some channels not used for ICA decomposition are used for rereferencing → the ICA decomposition has been removed
    % -> 이과정에서 EEG.icaweights은 지워짐
    append_report_row(log, cur_dname, cfg);
end
fprintf("Done\n")

Generate_perform_metrics

%% HELPERS %%
function log = init_log()
% INIT_LOG  MADE 리포트용 로그 struct를 기본값으로 초기화
%
%   log = init_log()
%
%   각 subject 처리 시작 시 호출한다. 모든 리포트 필드를 미리 기본값으로
%   채워두므로, subject가 중간에 skip(포기)되더라도 필드가 비지 않는다.
%   → 맨 끝에서 table을 만들거나 CSV로 append할 때 길이/필드 불일치가 안 생긴다.
%
%   필드 목록은 최종 report table의 열과 1:1로 대응한다.

    log.skip                                   = false;   % true면 이 subject 포기
    log.reference_used_for_faster              = 'NA';    % FASTER가 쓴 reference 라벨
    log.faster_bad_channels                    = '0';     % FASTER bad channel 목록
    log.ica_preparation_bad_channels           = '0';     % ICA 준비 단계 bad channel
    log.length_ica_data                        = 0;       % ICA 학습에 쓴 데이터 길이
    log.ica_data_rank                          = 0;       % ICA 데이터의 rank (=IC 수)
    log.total_ICs                              = 0;       % 총 IC 수
    log.ICs_removed                            = '0';     % 제거된 IC 목록
    log.adjust_performed                       = 'NA';    % ADJUST 실행 여부 (yes/skipped_rank/no_ica)
    log.total_epochs_before_artifact_rejection = 0;       % artifact 제거 전 epoch 수
    log.total_epochs_after_artifact_rejection  = 0;       % artifact 제거 후 epoch 수
    log.n_epoch_go_hit                         = 0;       % Go-Hit (CRN 대상)
    log.n_epoch_go_miss                        = 0;       % Go-miss
    log.n_epoch_nogo_fa                        = 0;       % NoGo-FA (ERN 대상)
    log.n_epoch_nogo_cr                        = 0;       % NoGo-CR
    log.n_epoch_uncat                          = 0;       % 미분류
    log.total_channels_interpolated            = 0;       % 보간된 채널 수
end


function append_report_row(log, data_name, cfg)
% APPEND_REPORT_ROW  한 subject의 log를 report CSV에 한 줄 append 하는 최소 스텁.
%
%   log struct의 각 필드를 한 행짜리 table로 만들어 output_location/MADE_report.csv 에
%   이어붙인다. 파일이 없으면 헤더 포함해서 새로 만든다.

    report_path = [cfg.output_location filesep cfg.report_name];

    % cell 배열 필드(reference_used_for_faster 등)는 문자열로 평탄화
    ref_used = log.reference_used_for_faster;
    if iscell(ref_used)
        ref_used = strjoin(ref_used, ',');
    end

    row = table( ...
        string(data_name), ...
        string(ref_used), ...
        string(log.faster_bad_channels), ...
        string(log.ica_preparation_bad_channels), ...
        log.length_ica_data, ...
        log.ica_data_rank, ...
        log.total_ICs, ...
        string(log.ICs_removed), ...
        string(log.adjust_performed), ...
        log.total_epochs_before_artifact_rejection, ...
        log.total_epochs_after_artifact_rejection, ...
        log.n_epoch_go_hit, ...
        log.n_epoch_go_miss, ...
        log.n_epoch_nogo_fa, ...
        log.n_epoch_nogo_cr, ...
        log.n_epoch_uncat, ...
        log.total_channels_interpolated, ...
        'VariableNames', { ...
            'filename', 'reference_used_for_faster', 'faster_bad_channels', ...
            'ica_preparation_bad_channels', 'length_ica_data', 'ica_data_rank', 'total_ICs', ...
            'ICs_removed', 'adjust_performed', 'total_epochs_before_artifact_rejection', ...
            'total_epochs_after_artifact_rejection', ...
            'n_epoch_go_hit', 'n_epoch_go_miss', 'n_epoch_nogo_fa', 'n_epoch_nogo_cr', 'n_epoch_uncat', ...
            'total_channels_interpolated'});

    if isfile(report_path)
        writetable(row, report_path, 'WriteMode', 'append');
    else
        writetable(row, report_path);
    end
end