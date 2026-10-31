% Wrapper MATLAB pour lancer Kolomvakis (BoolMF / CombineBoolMFs) sur
% une matrice au format Loom (CSV 0/1/-1) et emettre le CSV unifie.
%
% Necessite : MATLAB + Gurobi (avec licence academique).
% Sur Narval : module load matlab gurobi
%
% Usage (via matlab -batch) :
%   matlab -batch "run_kolomvakis('<input.csv>', <r>, <timeout_s>, ...
%                                 '<output.csv>', '<trace_out>', <seed>, '<regime>')"

function run_kolomvakis(input_csv, r, timeout_s, output_csv, trace_out, seed, regime)
    if nargin < 5, trace_out = ''; end
    if nargin < 6, seed = 42; end
    if nargin < 7, regime = 'raw'; end

    fprintf('[kolomvakis] input=%s r=%d timeout=%.1fs seed=%d\n', input_csv, r, timeout_s, seed);

    % Charger la matrice CSV (0/1/-1)
    M = csvread(input_csv);
    [m, n] = size(M);
    fprintf('[kolomvakis] M=%dx%d\n', m, n);

    % Ignorer les -1 en les mettant a 0 dans la matrice d'entree, mais
    % en conservant leur position pour l'evaluation d'erreur.
    mask_known = (M ~= -1);
    X = M; X(~mask_known) = 0;
    X = double(X);

    % Setup Kolomvakis (BooleanMF)
    kolomvakis_root = fullfile(fileparts(mfilename('fullpath')), ...
                               '..', '..', 'third_party', 'BooleanMF');
    addpath(genpath(kolomvakis_root));

    rng(seed);

    options = struct();
    options.totaltime = timeout_s;
    options.P = double(mask_known);  % ponderation : 0 pour les cellules -1
    options.contrRNG = 1;

    t0 = tic;
    trace_data = [];  % [time_ms err]
    try
        [W, H, err_wh, ~, ~] = CombineBoolMFs(X, r, options);
    catch ME
        fprintf('[kolomvakis] ERREUR : %s\n', ME.message);
        W = zeros(m, r); H = zeros(r, n);
        err_wh = sum(X(:) == 1);
    end
    total_ms = toc(t0) * 1000;

    % Calculer erreur reelle (ignorer les -1)
    prod = min(1, W * H);
    err_matrix = mask_known & (prod ~= (M == 1));
    final_err = sum(err_matrix(:));

    % Time-to-zero : approximation (CombineBoolMFs ne rapporte pas les
    % improvements intermediaires ; on met -1 si != 0, sinon total_ms).
    if final_err == 0
        time_to_zero_ms = total_ms;
    else
        time_to_zero_ms = -1;
    end

    fprintf('[kolomvakis] final_err=%d total_ms=%.1f\n', final_err, total_ms);

    % Trace : si dispo, une seule entree (init + final)
    if ~isempty(trace_out)
        init_err = sum(X(:) == 1);
        fid = fopen(trace_out, 'w');
        fprintf(fid, 'time_ms err\n');
        fprintf(fid, '0.0 %d\n', init_err);
        fprintf(fid, '%.3f %d\n', total_ms, final_err);
        fclose(fid);
    end

    % Emettre CSV au schema unifie
    [~, base, ~] = fileparts(input_csv);
    dataset = [base '.csv'];
    header = ['method,dataset,k,init_errors,phase1_errors,phase1_last_improve_ms,' ...
              'phase2_errors,phase2_last_improve_ms,final_errors,iterations,time_ms,' ...
              'wls_total_flips,wls_total_time_ms,wls_flips_per_sec,verif_errors,' ...
              'm,n,regime,tool,seed,time_to_zero_ms,trace_file'];
    init_err = sum(X(:) == 1);
    fid = fopen(output_csv, 'w');
    fprintf(fid, '%s\n', header);
    fprintf(fid, 'kolomvakis,%s,%d,%d,%d,%.1f,%d,%.1f,%d,0,%.1f,0,0.0,0.0,-1,%d,%d,%s,kolomvakis,%d,%.1f,%s\n', ...
        dataset, r, init_err, ...
        final_err, (time_to_zero_ms > 0) * time_to_zero_ms - (time_to_zero_ms <= 0), ...
        final_err, total_ms, ...
        final_err, total_ms, ...
        m, n, regime, seed, ...
        (time_to_zero_ms > 0) * time_to_zero_ms - (time_to_zero_ms <= 0), trace_out);
    fclose(fid);
    fprintf('[csv] %s\n', output_csv);
end
