function [EEG, log, state] = get_MADE_filtered_data(data_name, cfg, log, OVERLAP)
    
    state = struct();
    state.channels_analysed = [];
    state.ref_chan          = [];
    state.all_chan_bad_FAST = 0;

    output_dir = [cfg.output_location filesep 'filtered_data'];
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    [~, stem, ~] = fileparts(data_name);
    output_name = [output_dir filesep stem '_filtered_data' cfg.output_format];
    
    % 이미 filtering된 파일 있으면 로드하고 반환
    if ~OVERLAP && isfile(output_name)
        fprintf("\nPassing %s...(already filtered)\n", data_name)
        % Load saved filtered data file
        if strcmp(cfg.output_format, '.set')
            EEG = pop_loadset('filename', [stem '_filtered_data.set'], 'filepath', output_dir);
            if isfield(EEG, 'etc') && isfield(EEG.etc, 'MADE_log')
                log = EEG.etc.MADE_log; % log 복원
            else
                warning('MADE_log not found in %s. Using default log.', [stem '_filtered_data.set']);
            end
        elseif strcmp(cfg.output_format, '.mat')
            S = load([output_dir filesep stem '_filtered_data.mat'], 'EEG', 'log');
            EEG = S.EEG;
            if isfield(S, 'log')
                log = S.log; % log 복원       
            else
                warning('log not found in %s. Using default log.', [stem '_filtered_data.mat']);
            end
        end
        EEG = eeg_checkset(EEG);
        log.skip = false;
        return; 
    end

    fprintf("\nFiltering %s...\n", data_name)
    % STEP 1: Import
    EEG = mff_import([cfg.rawdata_location filesep data_name]);
    EEG = eeg_checkset(EEG);

    % STEP 1b: Type field 문자열화
    for atm = 1:length(EEG.event)
        if isnumeric(EEG.event(atm).type)
            EEG.event(atm).type = num2str(EEG.event(atm).type);
        end
    end

    % STEP 2: Channel locations
    EEG = pop_chanedit(EEG, 'lookup', cfg.channel_locations);
    EEG = eeg_checkset(EEG);
    if size(EEG.data, 1) ~= length(EEG.chanlocs)
        error('The size of the data does not match with channel numbers.');
    end

    % STEP 3: Time offset
    if cfg.adjust_time_offset == 1
        if cfg.filter_timeoffset ~= 0
            for aafto = 1:length(EEG.event)
                EEG.event(aafto).latency = EEG.event(aafto).latency ...
                    + (cfg.filter_timeoffset/1000)*EEG.srate;
            end
        end
        if cfg.stimulus_timeoffset ~= 0
            for sto = 1:length(EEG.event)
                for sm = 1:length(cfg.stimulus_markers)
                    if strcmp(EEG.event(sto).type, cfg.stimulus_markers{sm})
                        EEG.event(sto).latency = EEG.event(sto).latency ...
                            + (cfg.stimulus_timeoffset/1000)*EEG.srate;
                    end
                end
            end
        end
        if cfg.response_timeoffset ~= 0
            for rto = 1:length(EEG.event)
                for rm = 1:length(cfg.response_markers)
                    if strcmp(EEG.event(rto).type, cfg.response_markers{rm})
                        EEG.event(rto).latency = EEG.event(rto).latency ...
                            - (cfg.response_timeoffset/1000)*EEG.srate;
                    end
                end
            end
        end
    end

    % STEP 4: Resample
    if cfg.down_sample == 1
        if floor(cfg.sampling_rate) > EEG.srate
            error('Sampling rate cannot be higher than recorded sampling rate');
        elseif floor(cfg.sampling_rate) ~= EEG.srate
            EEG = pop_resample(EEG, cfg.sampling_rate);
            EEG = eeg_checkset(EEG);
        end
    end

    % STEP 5: 외곽채널 제거
    chans_labels = cell(1, EEG.nbchan);
    for i = 1:EEG.nbchan
        chans_labels{i} = EEG.chanlocs(i).labels;
    end
    if cfg.delete_outerlayer == 1
        [~, chansidx] = ismember(cfg.outerlayer_channel, chans_labels);
        outerlayer_channel_idx = chansidx(chansidx ~= 0);
        if isempty(outerlayer_channel_idx)
            error(['None of the outer layer channels present in channel locations of data.' ...
                ' Make sure outer layer channels are present in channel labels of data (EEG.chanlocs.labels).']);
        else
            EEG = pop_select(EEG, 'nochannel', outerlayer_channel_idx);
            EEG = eeg_checkset(EEG);
        end
    end

    % STEP 6: Filter
    high_transband = cfg.highpass;
    low_transband  = 10;
    hp_fl_order = 3.3 / (high_transband / EEG.srate);
    lp_fl_order = 3.3 / (low_transband  / EEG.srate);
    if mod(floor(hp_fl_order),2) == 0
        hp_fl_order = floor(hp_fl_order);
    else
        hp_fl_order = floor(hp_fl_order)+1;
    end
    if mod(floor(lp_fl_order),2) == 0
        lp_fl_order = floor(lp_fl_order)+2;
    else
        lp_fl_order = floor(lp_fl_order)+1;
    end
    high_cutoff = cfg.highpass/2;
    low_cutoff  = cfg.lowpass + (low_transband/2);

    EEG = eeg_checkset(EEG);
    EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass', ...
        'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
    EEG = eeg_checkset(EEG);
    EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass', ...
        'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
    EEG = eeg_checkset(EEG);

    % STEP 7: FASTER bad channels
    ref_chan = []; FASTbadChans = []; all_chan_bad_FAST = 0;
    ref_chan = find(any(EEG.data, 2)==0);

    if numel(ref_chan) > 1
        error(['There are more than 1 zeroed channel (i.e. zero value throughout recording) in data.' ...
            ' Only reference channel should be zeroed channel. Delete the zeroed channel/s which is not reference channel.']);

    elseif numel(ref_chan) == 1
        list_properties = channel_properties(EEG, 1:EEG.nbchan, ref_chan);
        FASTbadIdx  = min_z(list_properties);
        FASTbadChans = find(FASTbadIdx==1);
        FASTbadChans = FASTbadChans(FASTbadChans ~= ref_chan);
        log.reference_used_for_faster = {EEG.chanlocs(ref_chan).labels};
        EEG = eeg_checkset(EEG);
        state.channels_analysed = EEG.chanlocs;   % 보간용 전체 채널 위치 보관

    elseif numel(ref_chan) == 0
        warning('Reference channel is not present in data. Cz channel will be used as reference channel.');
        ref_chan = find(strcmp({EEG.chanlocs.labels}, 'Cz'));
        EEG_copy = EEG;
        EEG_copy = pop_reref(EEG_copy, ref_chan, 'keepref', 'on');
        EEG_copy = eeg_checkset(EEG_copy);
        list_properties = channel_properties(EEG_copy, 1:EEG_copy.nbchan, ref_chan);
        FASTbadIdx  = min_z(list_properties);
        FASTbadChans = find(FASTbadIdx==1);
        state.channels_analysed = EEG.chanlocs;
        log.reference_used_for_faster = {EEG.chanlocs(ref_chan).labels};
    end

    % 전멸 판정: 모든(또는 reference 제외 전부) 채널이 bad
    if numel(FASTbadChans)==EEG.nbchan || numel(FASTbadChans)+1==EEG.nbchan
        all_chan_bad_FAST = 1;
        warning(['No usable data for datafile ', data_name]);
        % all-bad 시점의 log를 완성해서 데이터와 함께 저장 (아래 skip 처리와 동일 값)
        log.faster_bad_channels                     = '0';
        log.ica_preparation_bad_channels            = '0';
        log.length_ica_data                         = 0;
        log.total_ICs                               = 0;
        log.ICs_removed                             = '0';
        log.total_epochs_before_artifact_rejection  = 0;
        log.total_epochs_after_artifact_rejection   = 0;
        log.total_channels_interpolated             = 0;
        log.skip                                    = true;
        if strcmp(cfg.output_format, '.set')
            EEG = eeg_checkset(EEG);
            EEG.etc.MADE_log = log;   % .set에는 EEG.etc에 실어서 저장
            EEG = pop_editset(EEG, 'setname', [stem '_no_usable_data_all_bad_channels']);
            EEG = pop_saveset(EEG, 'filename', [stem '_no_usable_data_all_bad_channels.set'], 'filepath', output_dir);
        elseif strcmp(cfg.output_format, '.mat')
            save([output_dir filesep stem '_no_usable_data_all_bad_channels.mat'], 'EEG', 'log');  % log 같이 저장
        end
    else
        % bad channel 제거
        EEG = pop_select(EEG, 'nochannel', FASTbadChans);
        EEG = eeg_checkset(EEG);
        if numel(ref_chan) == 1
            ref_chan = find(any(EEG.data, 2)==0);
            EEG = pop_select(EEG, 'nochannel', ref_chan);  % reference 채널 제거
        end
    end

    % faster_bad_channels 로그
    if numel(FASTbadChans) == 0
        log.faster_bad_channels = '0';
    else
        log.faster_bad_channels = num2str(FASTbadChans');
    end

    % state 갱신
    state.ref_chan          = ref_chan;
    state.all_chan_bad_FAST = all_chan_bad_FAST;
    state.FASTbadChans      = FASTbadChans;   % 뒤 STEP(채널 보간 수 집계)에서 사용

    % 전멸이면: 로그를 0으로 채우고 skip 플래그 세워서 반환 (continue 대체)
    if all_chan_bad_FAST == 1
        log.faster_bad_channels                     = '0';
        log.ica_preparation_bad_channels            = '0';
        log.length_ica_data                         = 0;
        log.total_ICs                               = 0;
        log.ICs_removed                             = '0';
        log.total_epochs_before_artifact_rejection  = 0;
        log.total_epochs_after_artifact_rejection   = 0;
        log.total_channels_interpolated             = 0;
        log.skip = true;
        return;   % 메인 루프에서 log.skip 보고 continue
    end

    % 저장 전에 log를 filtered 시점 상태로 확정 (로드 시 이 시점 그대로 복원됨)
    log.skip = false;

    % interim 저장 (옵션)
    if cfg.save_interim_result == 1
        % log를 데이터와 함께 저장
        if strcmp(cfg.output_format, '.set')
            EEG = eeg_checkset(EEG);
            EEG.etc.MADE_log = log;   % [B] .set에는 EEG.etc에 실어서 저장
            EEG = pop_editset(EEG, 'setname', [stem '_filtered_data']);
            EEG = pop_saveset(EEG, 'filename', [stem '_filtered_data.set'], 'filepath', output_dir);
        elseif strcmp(cfg.output_format, '.mat')
            save([output_dir filesep stem '_filtered_data.mat'], 'EEG', 'log');  % [B] log 같이 저장
        end
    end
end   