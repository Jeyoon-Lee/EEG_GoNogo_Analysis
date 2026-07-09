function EEG = Add_block_event(EEG, cur_dname, cfg)
% ADD_BLOCK_EVENT  단일 EEG의 EEG.event에 block/trial 번호 필드를 추가
%
%   EEG = Add_block_event(EEG)                 % 태깅만
%   EEG = Add_block_event(EEG, cur_dname, cfg) % 태깅 + filtered_data 파일에 저장
%
%   심는 필드 3개:
%     - block           : 몇 번째 block인지 (1부터)
%     - trial_in_block  : 그 block 안에서 몇 번째 trial (block마다 1로 리셋)
%     - trial_global    : 전체 통틀어 몇 번째 trial (1부터 쭉)
%
%   기준점:
%   trial 경계 : 'bgin' (모든 trial에 있으므로 이걸로 셈. bgin 하나 = trial 하나)
%   block 경계 : 직전 'TRSP' -> 다음 'bgin' gap이 gap_thresh_samp 초과 시 새 block.
%                (검증된 계산 방식: TRSP(k) -> bgin(k+1) 간격)
%
%   첫 bgin 이전 이벤트는 block=0, trial=0
%
%   계산이 불가능한 경우(이벤트/bgin 없음, bgin·TRSP 개수 불일치)에는 경고만 내고
%   EEG를 그대로(필드 미추가) 반환 (이 경우 저장도 하지 않음)
%
%   cur_dname, cfg가 주어지고 cfg.save_interim_result == 1이면, 태깅 성공 후
%   get_MADE_filtered_data가 저장한 것과 동일한 위치/이름의 filtered_data 파일을
%   덮어써서 block/trial 필드를 디스크에도 남김

    %% ===== PARAMETERS =====
    trial_marker    = 'bgin';
    end_marker      = 'TRSP';
    new_fields      = {'block', 'trial_in_block', 'trial_global'};
    gap_thresh_samp = 1000;   % block 경계 gap 임계(샘플). trial 간<10ms, block 간>14초. <- 직접 확인함
    expected_trials = 320;    % 8 block x 40 trial

    % 이미 field 존재할 경우 제거
    present = new_fields(isfield(EEG.event, new_fields));
    if ~isempty(present)
        EEG.event = rmfield(EEG.event, present);
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end

    nEv = numel(EEG.event);
    if nEv == 0
        warning('Add_block_event: 이벤트 없음. block/trial 필드 미추가.');
        return;
    end

    types = {EEG.event.type};
    lat   = [EEG.event.latency];

    bginIdx = find(strcmp(types, trial_marker));
    TRSPIdx = find(strcmp(types, end_marker));

    nBgin = numel(bginIdx);
    nTRSP = numel(TRSPIdx);

    %  gap 계산이 불가능한 경우 -> 필드 미추가 후 반환
    if nBgin == 0
        warning('Add_block_event: bgin이 0개. block/trial 필드 미추가.');
        return;
    end
    if nBgin ~= nTRSP
        warning('Add_block_event: bgin(%d)과 TRSP(%d) 개수 불일치. block/trial 필드 미추가.', nBgin, nTRSP);
        return;
    end

    % 320 trials이 아니어도 처리하되 flag만
    if nBgin ~= expected_trials
        warning('Add_block_event: bgin이 %d개 (기대 %d). 처리는 하되 확인 필요.', nBgin, expected_trials);
    end

    % 각 bgin(=trial)의 block 번호 계산 (검증된 벡터 방식)
    % TRSP(k) -> bgin(k+1) gap
    gaps = lat(bginIdx(2:end)) - lat(TRSPIdx(1:end-1));   % 길이 = nBgin-1
    is_new_block  = gaps > gap_thresh_samp;               % block 경계 여부
    block_of_bgin = [1, 1 + cumsum(is_new_block)];        % 각 trial의 block 번호 (길이 nBgin)

    % 각 trial의 block 내 trial 번호
    tinb_of_bgin = zeros(1, nBgin);
    for b = 1:max(block_of_bgin)
        mask = (block_of_bgin == b);
        tinb_of_bgin(mask) = 1:sum(mask);
    end
    tglob_of_bgin = 1:nBgin;   % 전체 통산 번호

    % 이벤트별로 소속 trial의 번호 물려주기
    % -> 각 이벤트는 "자기 앞(또는 자기)의 가장 가까운 bgin"이 속한 trial에 속함
    block_num = zeros(1, nEv);
    tinb_num  = zeros(1, nEv);
    tglob_num = zeros(1, nEv);

    cur = 0;   % 현재까지 지나온 bgin 순번 (0 = 아직 첫 bgin 전)
    for k = 1:nEv
        if strcmp(types{k}, trial_marker)
            cur = cur + 1;
        end
        if cur >= 1
            block_num(k) = block_of_bgin(cur);
            tinb_num(k)  = tinb_of_bgin(cur);
            tglob_num(k) = tglob_of_bgin(cur);
        end
        % cur==0 (첫 bgin 이전)이면 0 유지
    end

    % EEG.event에 필드 삽입
    for k = 1:nEv
        EEG.event(k).block          = block_num(k);
        EEG.event(k).trial_in_block = tinb_num(k);
        EEG.event(k).trial_global   = tglob_num(k);
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');

    % summary
    nBlocks = max(block_num);
    fprintf('  block 수: %d, 총 trial 수: %d\n', nBlocks, max(tglob_num));
    if nBlocks ~= 8
        warning('Add_block_event: block이 %d개 (기대 8).', nBlocks);
        for b = 1:nBlocks
            fprintf('    block %d: %d trials\n', b, sum(block_of_bgin==b));
        end
    end

    % filtered_data 파일에 태깅 결과 저장
    % get_MADE_filtered_data와 동일한 위치/이름으로 덮어써서 trial_global 등
    % 필드가 디스크 파일에도 남게 함
    if nargin >= 3 && isfield(cfg, 'save_interim_result') && cfg.save_interim_result == 1
        output_dir = [cfg.output_location filesep 'filtered_data'];
        [~, stem, ~] = fileparts(cur_dname);
        if strcmp(cfg.output_format, '.set')
            EEG = eeg_checkset(EEG);
            EEG = pop_saveset(EEG, 'filename', [stem '_filtered_data.set'], 'filepath', output_dir);
        elseif strcmp(cfg.output_format, '.mat')
            matpath = [output_dir filesep stem '_filtered_data.mat'];
            if isfile(matpath)
                S = load(matpath, 'log');   % 기존 log 보존해서 함께 저장
                log = S.log; %#ok<NASGU>
            else
                log = []; %#ok<NASGU>
            end
            save(matpath, 'EEG', 'log');
        end
    end
end
