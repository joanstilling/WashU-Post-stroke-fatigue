%% ===============================================================
% - FC (fs86, 86×86):
%     * Whole brain mean FC (upper triangle)
%     * Within-network FC (Yeo9)
%     * Between-network FC (Yeo9)
%     * ROI-of-interest connectivity:
%         - ROI -> network mean FC
%         - (optional) ROI -> ROI FC (86 targets)
%         - FC node strength
% - ALFF/fALFF
%     *ROI level 
%     *Network mean (Yeo9)
%
% - ChaCo (fs86 order):
%     * Robustly loads ANY 86×1 / 1×86 vector OR 86×86 matrix, regardless of variable name
%     * Exports:
%         - Network mean ChaCo (Yeo9)
%         - (optional) ROI-level ChaCo values
%
% Timepoint parsing:
%   sub-fcs002amc...   -> acute_arm_1
%   sub-fcs002amc2...  -> 3month_arm_1
%   (also supports older a/c/c2)
% ChaCo filenames like FCS_115_... assumed acute_arm_1
% ===============================================================

clear; clc;

%% ---------------- CONFIG ----------------
cfg.dataDir   = '/Users/...';

cfg.patternFC    = '*FCcorr*.mat';                 % FC files
cfg.patternChaCo = '*chaco*fs86*mean*.mat';        % adjust if needed (examples below)
cfg.patternALFF  = '*_alff_fs86_ts.mat';
cfg.patternFALFF = '*_falff_fs86_ts.mat';
cfg.export_roi_edge_fc = true;   % export all upper-triangle ROI<->ROI FC edges (3655/file)
cfg.roi_edge_use_indices_only = false;  % if true: ROI01<->ROI02 ; if false: ROI01_Name<->ROI02_Name

% Optional: export ROI-level ALFF/fALFF (can be big: 86 rows per file)
cfg.export_alff_roi_values  = true;
cfg.export_falff_roi_values = true;

% Always export network means (9 rows per file) unless you turn this off
cfg.export_alff_net_means  = true;
cfg.export_falff_net_means = true;

cfg.outFile = fullfile(cfg.dataDir, 'regionFC_and_ChaCo_alff_falff_features_long.csv');


% DLPFC ROI(s) of interest (fs86 ROI names)
cfg.dlpfc_rois = [ ...
    "ctx-lh-rostralmiddlefrontal", ...
    "ctx-rh-rostralmiddlefrontal" ...
];
cfg.export_seed_pair_edges = true;   % turn on seed<->seed summary edges

% ---- Secondary seeds (bilateral) ----
cfg.secondary_seed_sets = struct();
cfg.secondary_seed_sets.THalamus  = ["Left-Thalamus-Proper","Right-Thalamus-Proper"];
cfg.secondary_seed_sets.Precuneus = ["ctx-lh-precuneus","ctx-rh-precuneus"];
cfg.secondary_seed_sets.Insula    = ["ctx-lh-insula","ctx-rh-insula"];
cfg.secondary_seed_sets.CaudalACC = ["ctx-lh-caudalanteriorcingulate","ctx-rh-caudalanteriorcingulate"]; % optional

% Optional exports (can get big)
cfg.export_secondary_seed_to_all_rois = false;  % set true if you want seed->ROI for all targets


% Optional exports (can get big)
cfg.export_dlpfc_to_all_rois = true;   % DLPFC->ROI (85 values per file)
cfg.export_chaco_roi_values  = false;  % ROI ChaCo (86 values per file). Keep false unless needed.

%% ---------------- fs86 ROI names ----------------
region_names = [ ...
"Left-Cerebellum-Cortex", "Left-Thalamus-Proper", "Left-Caudate", ...
"Left-Putamen", "Left-Pallidum", "Left-Hippocampus", "Left-Amygdala", ...
"Left-Accumbens-area", "Left-VentralDC", "Right-Cerebellum-Cortex", ...
"Right-Thalamus-Proper", "Right-Caudate", "Right-Putamen", ...
"Right-Pallidum", "Right-Hippocampus", "Right-Amygdala", ...
"Right-Accumbens-area", "Right-VentralDC", "ctx-lh-bankssts", ...
"ctx-lh-caudalanteriorcingulate", "ctx-lh-caudalmiddlefrontal", ...
"ctx-lh-cuneus", "ctx-lh-entorhinal", "ctx-lh-fusiform", ...
"ctx-lh-inferiorparietal", "ctx-lh-inferiortemporal", ...
"ctx-lh-isthmuscingulate", "ctx-lh-lateraloccipital", ...
"ctx-lh-lateralorbitofrontal", "ctx-lh-lingual", ...
"ctx-lh-medialorbitofrontal", "ctx-lh-middletemporal", ...
"ctx-lh-parahippocampal", "ctx-lh-paracentral", ...
"ctx-lh-parsopercularis", "ctx-lh-parsorbitalis", ...
"ctx-lh-parstriangularis", "ctx-lh-pericalcarine", ...
"ctx-lh-postcentral", "ctx-lh-posteriorcingulate", ...
"ctx-lh-precentral", "ctx-lh-precuneus", ...
"ctx-lh-rostralanteriorcingulate", "ctx-lh-rostralmiddlefrontal", ...
"ctx-lh-superiorfrontal", "ctx-lh-superiorparietal", ...
"ctx-lh-superiortemporal", "ctx-lh-supramarginal", "ctx-lh-frontalpole", ...
"ctx-lh-temporalpole", "ctx-lh-transversetemporal", "ctx-lh-insula", ...
"ctx-rh-bankssts", "ctx-rh-caudalanteriorcingulate", ...
"ctx-rh-caudalmiddlefrontal", "ctx-rh-cuneus", "ctx-rh-entorhinal", ...
"ctx-rh-fusiform", "ctx-rh-inferiorparietal", ...
"ctx-rh-inferiortemporal", "ctx-rh-isthmuscingulate", ...
"ctx-rh-lateraloccipital", "ctx-rh-lateralorbitofrontal", ...
"ctx-rh-lingual", "ctx-rh-medialorbitofrontal", ...
"ctx-rh-middletemporal", "ctx-rh-parahippocampal", "ctx-rh-paracentral", ...
"ctx-rh-parsopercularis", "ctx-rh-parsorbitalis", ...
"ctx-rh-parstriangularis", "ctx-rh-pericalcarine", "ctx-rh-postcentral", ...
"ctx-rh-posteriorcingulate", "ctx-rh-precentral", "ctx-rh-precuneus", ...
"ctx-rh-rostralanteriorcingulate", "ctx-rh-rostralmiddlefrontal", ...
"ctx-rh-superiorfrontal", "ctx-rh-superiorparietal", ...
"ctx-rh-superiortemporal", "ctx-rh-supramarginal", "ctx-rh-frontalpole", ...
"ctx-rh-temporalpole", "ctx-rh-transversetemporal", "ctx-rh-insula" ...
]';

nROI = numel(region_names);

%% ---------------- Yeo9 labels (1..9) ----------------
network_labels = [ ...
9,8,8,8,8,8,8,8,8,9,8,8,8,8,8,8,8,8, ...
7,4,7,1,5,1,7,5,7,1,5,1,5,7,1,2,7,7,7, ...
1,2,4,2,7,7,6,7,3,2,4,5,5,2,4,2,4,6,1,5, ...
1,7,5,7,1,5,1,5,7,1,2,6,7,6,1,2,4,2,7,7, ...
6,7,3,2,4,5,5,2,4 ...
]';

network_names = ["VIS","SOM","DAN","VAN","LIM","FP","DMN","SUB","CER"];
nNet = numel(network_names);

roiIdxByNet = cell(nNet,1);
for k = 1:nNet
    roiIdxByNet{k} = find(network_labels == k);
end

%% ---------------- DLPFC indices ----------------
dlpfc_idx = find(ismember(region_names, cfg.dlpfc_rois));
if isempty(dlpfc_idx)
    error('None of the requested DLPFC ROI(s) were found in region_names.');
end
fprintf('DLPFC ROIs:\n');
for i = 1:numel(dlpfc_idx)
    fprintf('  %s (index %d)\n', region_names(dlpfc_idx(i)), dlpfc_idx(i));
end

%% ---------------- Secondary seed indices ----------------
seed_names = fieldnames(cfg.secondary_seed_sets);
seed_idx_map = struct();

fprintf('\nSecondary seeds:\n');
for s = 1:numel(seed_names)
    nm = seed_names{s};
    rois = cfg.secondary_seed_sets.(nm);
    idx = find(ismember(region_names, rois));

    if isempty(idx)
        warning('Secondary seed %s: none of ROIs found: %s', nm, strjoin(cellstr(rois), ', '));
    else
        seed_idx_map.(nm) = idx;
        fprintf('  %s:\n', nm);
        for j = 1:numel(idx)
            fprintf('    %s (index %d)\n', region_names(idx(j)), idx(j));
        end
    end
end

%% ---------------- Unified seed index map (includes DLPFC) ----------------
all_seed_idx_map = seed_idx_map;  % secondary seeds you found

% add DLPFC as a named seed (so we can compute DLPFC<->others uniformly)
all_seed_idx_map.DLPFC = dlpfc_idx;

all_seed_names = fieldnames(all_seed_idx_map);

fprintf('\nAll seeds used for seed-pair edges:\n');
for s = 1:numel(all_seed_names)
    nm = all_seed_names{s};
    fprintf('  %s: %d nodes\n', nm, numel(all_seed_idx_map.(nm)));
end

%% ---------------- Precompute masks & names ----------------
maskROI = triu(true(nROI),1);

maskNet = triu(true(nNet),1);
[ni, nj] = find(maskNet);
netEdgeNames = strings(numel(ni),1);
for e = 1:numel(ni)
    netEdgeNames(e) = network_names(ni(e)) + "-" + network_names(nj(e));
end

%% ---------------- Find files (recursive) ----------------
[fcFiles, idsFC, eventsFC] = list_files_recursive(cfg.dataDir, cfg.patternFC);
[chFiles, idsCH, eventsCH] = list_files_recursive(cfg.dataDir, cfg.patternChaCo);
[alffFiles,  idsALFF,  eventsALFF]  = list_files_recursive(cfg.dataDir, cfg.patternALFF);
[falffFiles, idsFALFF, eventsFALFF] = list_files_recursive(cfg.dataDir, cfg.patternFALFF);

fprintf('Found %d ALFF files, %d fALFF files.\n', numel(alffFiles), numel(falffFiles));

fprintf('\nFound %d FC files, %d ChaCo files.\n', numel(fcFiles), numel(chFiles));

if ~isempty(fcFiles)
    fprintf('FC timepoint counts:\n');
    u = unique(eventsFC);
    for k = 1:numel(u)
        fprintf('  %s: %d\n', u(k), sum(eventsFC==u(k)));
    end
end

if ~isempty(chFiles)
    fprintf('ChaCo timepoint counts:\n');
    u = unique(eventsCH);
    for k = 1:numel(u)
        fprintf('  %s: %d\n', u(k), sum(eventsCH==u(k)));
    end
end

if ~isempty(alffFiles)
    fprintf('ALFF timepoint counts:\n');
    u = unique(eventsALFF);
    for k = 1:numel(u), fprintf('  %s: %d\n', u(k), sum(eventsALFF==u(k))); end
end

if ~isempty(falffFiles)
    fprintf('fALFF timepoint counts:\n');
    u = unique(eventsFALFF);
    for k = 1:numel(u), fprintf('  %s: %d\n', u(k), sum(eventsFALFF==u(k))); end
end
%% ---------------- Extract features ----------------
T_all = table();

%% ===== FC FEATURES =====
for s = 1:numel(fcFiles)
    fcFile = fullfile(fcFiles(s).folder, fcFiles(s).name);
    study_id   = idsFC(s);
    event_name = eventsFC(s);

    S = load(fcFile);

    % Robust FC matrix extraction (86×86), regardless of variable name
    M = find_numeric_matrix_in_struct(S, [nROI nROI], {'FCcorr','conmat','FC','R','Z','corrmat'});
    if isempty(M)
        warning('Skipping %s: no %dx%d FC matrix found.', fcFiles(s).name, nROI, nROI);
        continue;
    end

    M = double(M);
    M = (M + M.')/2;
    M(1:nROI+1:end) = 0;

    rows = table();

     % If it looks like correlation (mostly within [-1,1]), convert to Fisher-z
    if max(abs(M(:)), [], 'omitnan') <= 1.001
        M = atanh(M);
    end


    % A) Whole-brain mean FC
    wb_mean = mean(M(maskROI), 'omitnan');
    rows = [rows; make_row(study_id, event_name, "WB_FC", "wholebrain_mean", wb_mean)];

    % B) Within-network FC
    for k = 1:nNet
        idx = roiIdxByNet{k};
        if numel(idx) < 2, continue; end
        subM = M(idx, idx);
        subMask = triu(true(numel(idx)), 1);
        val = mean(subM(subMask), 'omitnan');
        rows = [rows; make_row(study_id, event_name, "NET_WITHIN_FC", "within_" + network_names(k), val)];
    end

    % C) Between-network FC
    NetMat = nan(nNet);
    for a = 1:nNet
        Ia = roiIdxByNet{a};
        if isempty(Ia), continue; end
        for b = a+1:nNet
            Ib = roiIdxByNet{b};
            if isempty(Ib), continue; end
            val = mean(M(Ia, Ib), 'all', 'omitnan');
            NetMat(a,b) = val;
            NetMat(b,a) = val;
        end
    end
    between_vals = NetMat(maskNet);
    for e = 1:numel(between_vals)
        rows = [rows; make_row(study_id, event_name, "NET_BETWEEN_FC", netEdgeNames(e), between_vals(e))];
    end

    % D) DLPFC -> network mean FC
    for k = 1:nNet
        idx = roiIdxByNet{k};
        val = mean(M(dlpfc_idx, idx), 'all', 'omitnan');
        rows = [rows; make_row(study_id, event_name, "DLPFC_TO_NET_FC", "DLPFC->" + network_names(k), val)];
    end

    % E) Optional: DLPFC -> ROI (all targets)
    if cfg.export_dlpfc_to_all_rois
        for r = 1:nROI
            if ismember(r, dlpfc_idx), continue; end
            val = mean(M(dlpfc_idx, r), 'all', 'omitnan');
            rows = [rows; make_row(study_id, event_name, "DLPFC_TO_ROI_FC", "DLPFC->" + region_names(r), val)];
        end
    end
    % ==========================================================
    % F) Secondary seed connectivity (seed -> network, optional seed -> ROI)
    % ==========================================================
    for s2 = 1:numel(seed_names)
        seed_nm = seed_names{s2};

        if ~isfield(seed_idx_map, seed_nm)
            continue; % seed not found in atlas
        end

        seed_idx = seed_idx_map.(seed_nm);

        % ---- Seed -> network mean FC ----
        for k = 1:nNet
            idx = roiIdxByNet{k};
            val = mean(M(seed_idx, idx), 'all', 'omitnan');
            rows = [rows; make_row(study_id, event_name, ...
                "SEED_TO_NET_FC", "SEED:" + seed_nm + "->" + network_names(k), val)];
        end

        % ---- Optional: Seed -> ROI FC (all targets) ----
        if cfg.export_secondary_seed_to_all_rois
            for r = 1:nROI
                if ismember(r, seed_idx), continue; end
                val = mean(M(seed_idx, r), 'all', 'omitnan');
                rows = [rows; make_row(study_id, event_name, ...
                    "SEED_TO_ROI_FC", "SEED:" + seed_nm + "->" + region_names(r), val)];
            end
        end
    end
    % ==========================================================
    % G) Seed <-> Seed summary edges (small set, very interpretable)
    % ==========================================================
    if isfield(cfg, 'export_seed_pair_edges') && cfg.export_seed_pair_edges

        % pairwise combinations of all seeds
        for a = 1:numel(all_seed_names)
            Aname = string(all_seed_names{a});
            Aidx  = all_seed_idx_map.(all_seed_names{a});

            for b = a+1:numel(all_seed_names)
                Bname = string(all_seed_names{b});
                Bidx  = all_seed_idx_map.(all_seed_names{b});

                % mean FC across all A x B edges (handles bilateral/multi-node seeds)
                val = mean(M(Aidx, Bidx), 'all', 'omitnan');

                rows = [rows; make_row(study_id, event_name, ...
                    "SEED_EDGE_FC", Aname + "<->" + Bname, val)];
            end
        end
    end
        % ----------------------------------------------
        % H) FC NODE STRENGTH (parcel-wise)
        % ----------------------------------------------
        Mns = M;
        Mns(1:nROI+1:end) = NaN;                 % exclude diagonal
        node_strength = mean(Mns, 2, 'omitnan');  % 86x1 (mean connectivity per node)
        
        for r = 1:nROI
            rows = [rows; make_row(study_id, event_name, ...
                "FC_NODE_STRENGTH", "NS_" + region_names(r), node_strength(r))];
        end 

    T_all = [T_all; rows]; %#ok<AGROW>
end
key = strcat(T_all.study_id,"|",T_all.event_name,"|",T_all.feature_type,"|",T_all.feature_name);
[ukey,~,ic] = unique(key);
counts = accumarray(ic, 1);
dupKeys = ukey(counts > 1);

fprintf("Duplicate feature rows (same id/event/type/name): %d\n", numel(dupKeys));
if ~isempty(dupKeys)
    disp(dupKeys(1:min(20,end)));
end

%% ===== ALFF FEATURES (fs86 ROI vector) =====
for s = 1:numel(alffFiles)
    f = fullfile(alffFiles(s).folder, alffFiles(s).name);
    study_id   = idsALFF(s);
    event_name = eventsALFF(s);

    S = load(f);

    % Robustly extract ROI vector (86x1) even if stored as 1x86, 86xT, Tx86
    v = find_roi_vector_in_struct(S, nROI, {'alff','ALFF','fs86_alff','alff_fs86','ts','data','X'});
    if isempty(v)
        warning('Skipping %s: no fs86-length ALFF vector found.', alffFiles(s).name);
        continue;
    end
    v = double(v(:));  % force column

    rows = table();

    % A) ROI-level ALFF
    if isfield(cfg,'export_alff_roi_values') && cfg.export_alff_roi_values
        for r = 1:nROI
            rows = [rows; make_row(study_id, event_name, "ALFF_ROI", "ALFF_" + region_names(r), v(r))];
        end
    end

    % B) Network mean ALFF (Yeo9)
    if ~isfield(cfg,'export_alff_net_means') || cfg.export_alff_net_means
        for k = 1:nNet
            idx = roiIdxByNet{k};
            val = mean(v(idx), 'omitnan');
            rows = [rows; make_row(study_id, event_name, "ALFF_NET_MEAN", "ALFF_" + network_names(k), val)];
        end
    end

    T_all = [T_all; rows]; %#ok<AGROW>
end

%% ===== fALFF FEATURES (fs86 ROI vector) =====
for s = 1:numel(falffFiles)
    f = fullfile(falffFiles(s).folder, falffFiles(s).name);
    study_id   = idsFALFF(s);
    event_name = eventsFALFF(s);

    S = load(f);

    v = find_roi_vector_in_struct(S, nROI, {'falff','fALFF','FALFF','fs86_falff','falff_fs86','ts','data','X'});
    if isempty(v)
        warning('Skipping %s: no fs86-length fALFF vector found.', falffFiles(s).name);
        continue;
    end
    v = double(v(:));

    rows = table();

    % A) ROI-level fALFF
    if isfield(cfg,'export_falff_roi_values') && cfg.export_falff_roi_values
        for r = 1:nROI
            rows = [rows; make_row(study_id, event_name, "FALFF_ROI", "fALFF_" + region_names(r), v(r))];
        end
    end

    % B) Network mean fALFF (Yeo9)
    if ~isfield(cfg,'export_falff_net_means') || cfg.export_falff_net_means
        for k = 1:nNet
            idx = roiIdxByNet{k};
            val = mean(v(idx), 'omitnan');
            rows = [rows; make_row(study_id, event_name, "FALFF_NET_MEAN", "fALFF_" + network_names(k), val)];
        end
    end

    T_all = [T_all; rows]; %#ok<AGROW>
end
%% ===== ChaCo FEATURES (EDGE MATRIX) =====
for s = 1:numel(chFiles)
    chFile = fullfile(chFiles(s).folder, chFiles(s).name);
    study_id   = idsCH(s);
    event_name = eventsCH(s);  % your parser sets acute for FCS_### by default

    S = load(chFile);
    % inside ChaCo loop, after parsing id/event:
    if startsWith(study_id, "amc_")
        continue; % skip controls for ChaCo export
    end

    % Expect NeMo chacoconn output: variable "C" (86x86)
    if isfield(S, 'C') && isnumeric(S.C) && all(size(S.C) == [nROI nROI])
        C = S.C;
    else
        % Fallback: search any 86x86 numeric matrix
        C = find_numeric_matrix_in_struct(S, [nROI nROI], {'C','chacoconn','ChaCoConn','chaco_conn','ChaCoMat','chaco'});
        if isempty(C)
            fprintf('Skipping %s: no 86x86 ChaCoConn matrix found.\n', chFiles(s).name);
            continue;
        end
    end

    % ---- CRITICAL: ensure full double (not sparse) ----
    C = full(double(C));

    % ---- Symmetry check (scalar) ----
    asym = max(abs(C - C.'), [], 'all');  % scalar

    if asym > 1e-6
        fprintf('NOTE: ChaCoConn not symmetric in %s (max asym=%g). Symmetrizing.\n', ...
            char(chFiles(s).name), asym);
        C = (C + C.') / 2;
    end

    % Diagonal is not meaningful for edge disruption; exclude it
    C(1:nROI+1:end) = NaN;

    rows = table();

    % =========================
    % A) ROI "burden" (node strength-like): mean disruption of edges incident to ROI
    % =========================
    roi_burden = mean(C, 2, 'omitnan');  % 86x1

    % Export ROI-level burden if desired
    if cfg.export_chaco_roi_values
        for r = 1:nROI
            rows = [rows; make_row(study_id, event_name, "CHACO_ROI_BURDEN", "ChaCoBurden_" + region_names(r), roi_burden(r))];
        end
    end

    % Network mean burden (network-level node disruption)
    for k = 1:nNet
        idx = roiIdxByNet{k};
        val = mean(roi_burden(idx), 'omitnan');
        rows = [rows; make_row(study_id, event_name, "CHACO_NET_MEAN", "ChaCoBurden_" + network_names(k), val)];
    end

    % =========================
    % B) Within-network mean edge disruption
    % =========================
    for k = 1:nNet
        idx = roiIdxByNet{k};
        if numel(idx) < 2, continue; end
        subM = C(idx, idx);
        subMask = triu(true(numel(idx)), 1);
        val = mean(subM(subMask), 'omitnan');
        rows = [rows; make_row(study_id, event_name, "CHACO_NET_WITHIN_MEAN", "within_" + network_names(k), val)];
    end

    % =========================
    % C) Between-network mean edge disruption
    % =========================
    NetMat = nan(nNet);
    for a = 1:nNet
        Ia = roiIdxByNet{a};
        if isempty(Ia), continue; end
        for b = a+1:nNet
            Ib = roiIdxByNet{b};
            if isempty(Ib), continue; end
            val = mean(C(Ia, Ib), 'all', 'omitnan');
            NetMat(a,b) = val;
            NetMat(b,a) = val;
        end
    end

    between_vals = NetMat(maskNet);
    for e = 1:numel(between_vals)
        rows = [rows; make_row(study_id, event_name, "CHACO_NET_BETWEEN_MEAN", "between_" + netEdgeNames(e), between_vals(e))];
    end

    T_all = [T_all; rows]; %#ok<AGROW>
end

% ---- DEBUG: confirm ChaCo feature types produced ----
isChaco = startsWith(string(T_all.feature_type), "CHACO");
fprintf('\nChaCo feature_type counts in T_all (before write):\n');
disp(groupsummary(T_all(isChaco,:), "feature_type"));

% Explicit counts (more robust than groupsummary in older MATLAB)
fprintf('CHACO_NET_WITHIN_MEAN rows: %d\n', sum(string(T_all.feature_type)=="CHACO_NET_WITHIN_MEAN"));
fprintf('CHACO_NET_BETWEEN_MEAN rows: %d\n', sum(string(T_all.feature_type)=="CHACO_NET_BETWEEN_MEAN"));
fprintf('CHACO_NET_MEAN rows: %d\n', sum(string(T_all.feature_type)=="CHACO_NET_MEAN"));


%% ---------------- Write output ----------------
writetable(T_all, cfg.outFile);
fprintf('\nWrote %s with %d rows.\n', cfg.outFile, height(T_all));

%% ===============================================================
% Row helper
function row = make_row(study_id, event_name, feature_type, feature_name, value)
    row = table(string(study_id), string(event_name), string(feature_type), string(feature_name), double(value), ...
        'VariableNames', {'study_id','event_name','feature_type','feature_name','value'});
end

%% ===============================================================
% Recursive file listing + ID/event parsing
function [fileList, ids, events] = list_files_recursive(dirPath, pattern)
    fileList = dir(fullfile(dirPath, '**', pattern));
    nF = numel(fileList);
    ids = strings(nF,1);
    events = strings(nF,1);
    for i = 1:nF
        [ids(i), events(i)] = extract_id_and_event(fileList(i).name);
    end
end

%% ===============================================================
% ID + event parser (FC "sub-..." + ChaCo "FCS_###...")
function [id, event_name] = extract_id_and_event(fname)
% Robust subject + event parser for WashU FC + ChaCo filenames.
%
% FC examples:
%   sub-fcs002amc_...  -> id=amc_002_fcs, event=visit_1_arm_2
%   sub-fcs002amc2_... -> id=amc_002_fcs, event=visit_2_arm_2
%   sub-fcs024a_...    -> id=fcs_024,     event=acute_arm_1
%   sub-fcs024c_...    -> id=fcs_024,     event=3month_arm_1
%   sub-fcs024c2_...   -> id=fcs_024,     event=1year_arm_1
%
% ChaCo example:
%   FCS_115_nemo_output_... -> id=fcs_115, event=acute_arm_1 (assumed)

    fname = string(fname);
    id = "";
    event_name = "";

    % ---------- Case 1: FC files with sub- prefix ----------
    tok = regexp(fname, '^sub-([^_]+)', 'tokens', 'once');
    if ~isempty(tok)
        raw = string(tok{1});   % e.g., fcs002amc2, fcs024a, fcs024c2

        m = regexp(raw, '^fcs(\d+)(amc2|amc|c2|a|c)?$', 'tokens', 'once');
        if isempty(m)
            warning('Could not parse subject token from filename: %s', fname);
            id = raw;
            event_name = "";
            return;
        end

        numStr = string(m{1});            % e.g., "002" or "24"
        numStr = pad(numStr, 3, 'left', '0');  % -> "024"

        suffix = "";
        if numel(m) >= 2 && ~isempty(m{2})
            suffix = string(m{2});        % amc, amc2, a, c, c2
        end

        % ID mapping
        if suffix == "amc" || suffix == "amc2"
            id = "amc_" + numStr + "_fcs";   % controls in REDCap: AMC_###_FCS
        else
            id = "fcs_" + numStr;            % stroke in REDCap: FCS_###
        end

        % Event mapping
        switch suffix
            case {"amc"}
                event_name = "visit_1_arm_2";
            case "amc2"
                event_name = "visit_2_arm_2";
            case "a"
                event_name = "acute_arm_1";
            case "c"
                event_name = "3month_arm_1";
            case "c2"
                event_name = "1year_arm_1";
          
        end

        return;
    end

    % ---------- Case 2: ChaCo files like FCS_115_... or AMC_001_FCS_... ----------
    m2 = regexp(fname, '^(AMC|FCS)_(\d+)', 'tokens', 'once');
    if ~isempty(m2)
        prefix = string(m2{1});  % AMC or FCS
        numStr = string(m2{2});
        numStr = pad(numStr, 3, 'left', '0');

        if prefix == "AMC"
            id = "amc_" + numStr + "_fcs";
        else
            id = "fcs_" + numStr;
        end

        % You told me ChaCo is acute-only
        event_name = "acute_arm_1";
        return;
    end

    % ---------- Fallback ----------
    warning('Could not parse ID/event from filename: %s', fname);
    id = "";
    event_name = "";
end 


%% ===============================================================
% Robust extractor: find numeric matrix of specified size in a loaded struct
function M = find_numeric_matrix_in_struct(S, targetSize, preferredNames)
    M = [];

    % 1) Try preferred names first
    for i = 1:numel(preferredNames)
        nm = preferredNames{i};
        if isfield(S, nm)
            X = S.(nm);
            if isnumeric(X) && ismatrix(X) && all(size(X) == targetSize)
                M = X;
                return;
            end
        end
    end

    % 2) Otherwise scan all fields for a matching numeric matrix
    fn = fieldnames(S);
    for k = 1:numel(fn)
        X = S.(fn{k});
        if isnumeric(X) && ismatrix(X) && all(size(X) == targetSize)
            M = X;
            return;
        end
    end
end

%% ===============================================================
% Robust extractor: find numeric vector of specified length in a loaded struct
function v = find_numeric_vector_in_struct(S, targetLen, preferredNames)
    v = [];

    % 1) Try preferred names first
    for i = 1:numel(preferredNames)
        nm = preferredNames{i};
        if isfield(S, nm)
            X = S.(nm);
            if isnumeric(X) && isvector(X) && numel(X) == targetLen
                v = X;
                return;
            end
        end
    end

    % 2) Otherwise scan all fields for any numeric vector of the right length
    fn = fieldnames(S);
    for k = 1:numel(fn)
        X = S.(fn{k});
        if isnumeric(X) && isvector(X) && numel(X) == targetLen
            v = X;
            return;
        end
    end
end
%% ===============================================================
% Robust extractor: find fs86 ROI vector even if stored as 86x1, 1x86, 86xT, or Tx86
function v = find_roi_vector_in_struct(S, nROI, preferredNames)
    v = [];

    % ---- 1) Preferred names first ----
    for i = 1:numel(preferredNames)
        nm = preferredNames{i};
        if isfield(S, nm)
            X = S.(nm);
            v = coerce_to_roi_vector(X, nROI);
            if ~isempty(v), return; end
        end
    end

    % ---- 2) Scan all fields ----
    fn = fieldnames(S);
    for k = 1:numel(fn)
        X = S.(fn{k});
        v = coerce_to_roi_vector(X, nROI);
        if ~isempty(v), return; end
    end
end

function v = coerce_to_roi_vector(X, nROI)
    v = [];
    if ~isnumeric(X) || isempty(X)
        return;
    end

    % Case A: already a vector of length nROI
    if isvector(X) && numel(X) == nROI
        v = double(X(:));
        return;
    end

    % Case B: matrix with one dimension = nROI (86xT or Tx86)
    if ismatrix(X)
        sz = size(X);
        if any(sz == nROI)
            % If 86xT -> mean across T (dim 2)
            if sz(1) == nROI
                v = mean(double(X), 2, 'omitnan');
                return;
            end
            % If Tx86 -> mean across T (dim 1)
            if sz(2) == nROI
                v = mean(double(X), 1, 'omitnan').';
                return;
            end
        end
    end
end