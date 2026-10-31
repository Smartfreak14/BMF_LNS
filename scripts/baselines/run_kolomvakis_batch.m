function run_kolomvakis_batch(datasets_csv, data_dir, output_dir, nseeds, timeout_override, regime)
% Batch : lance Kolomvakis (BoolMF+CombineBoolMFs) sur tous les datasets
% listes dans datasets_csv, NSEEDS graines chacun.
%
% Args :
%   datasets_csv     : chemin vers un CSV 'filename,k,timeout' (comme
%                      cluster/datasets_real.csv ou datasets_table2.csv)
%   data_dir         : dossier des matrices CSV (ex '.../data/Real')
%   output_dir       : dossier de sortie (cree si absent)
%   nseeds           : graines (defaut 5)
%   timeout_override : timeout uniforme secondes ; vide=[] utilise csv
%   regime           : etiquette 'raw'/'simpluni'/'simplexi' (defaut 'raw')
%
% Exemple :
%   run_kolomvakis_batch('/home/kiko/.../cluster/datasets_real.csv', ...
%                        '/home/kiko/.../data/Real', ...
%                        '/home/kiko/.../cluster/results_kolomvakis_real', ...
%                        5, [], 'raw')

if nargin < 4, nseeds = 5; end
if nargin < 5, timeout_override = []; end
if nargin < 6, regime = 'raw'; end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Charger les modules kolomvakis
here = fileparts(mfilename('fullpath'));
kolomvakis_root = fullfile(here, '..', '..', 'third_party', 'BooleanMF');
addpath(genpath(kolomvakis_root));
addpath(here);   % pour run_kolomvakis

% Parser datasets_csv (skip entete + lignes commentees)
fid = fopen(datasets_csv, 'r');
datasets = {};
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || line(1) == '#' || strcmp(line, 'filename,k,timeout_per_option_seconds')
        continue
    end
    parts = strsplit(line, ',');
    if length(parts) < 2, continue; end
    filename = strtrim(parts{1});
    k = str2double(strtrim(parts{2}));
    if length(parts) >= 3
        tout = str2double(strtrim(parts{3}));
    else
        tout = 300;
    end
    if ~isempty(timeout_override), tout = timeout_override; end
    datasets{end+1, 1} = {filename, k, tout};
end
fclose(fid);

fprintf('[batch] %d datasets x %d seeds = %d runs\n', ...
    length(datasets), nseeds, length(datasets)*nseeds);
fprintf('[batch] data_dir=%s\n', data_dir);
fprintf('[batch] output_dir=%s\n', output_dir);

n_ok = 0;
n_fail = 0;
for idx = 1:length(datasets)
    entry = datasets{idx};
    filename = entry{1};
    k = entry{2};
    tout = entry{3};

    input_path = fullfile(data_dir, filename);
    if ~exist(input_path, 'file')
        fprintf('[batch] [%d/%d] SKIP %s : absent\n', idx, length(datasets), filename);
        n_fail = n_fail + 1;
        continue
    end
    [~, base, ~] = fileparts(filename);
    fprintf('[batch] [%d/%d] %s k=%d timeout=%.1fs\n', ...
        idx, length(datasets), filename, k, tout);
    for seed = 1:nseeds
        out_csv = fullfile(output_dir, sprintf('%s_kolomvakis_s%d.csv', base, seed));
        trace_out = fullfile(output_dir, sprintf('%s_kolomvakis_s%d.trace', base, seed));
        try
            run_kolomvakis(input_path, k, tout, out_csv, trace_out, seed, regime);
            n_ok = n_ok + 1;
        catch ME
            fprintf('[batch]   [seed=%d] ERREUR : %s\n', seed, ME.message);
            n_fail = n_fail + 1;
        end
    end
end

fprintf('[batch] Fini : %d OK, %d echec/skip\n', n_ok, n_fail);
end
