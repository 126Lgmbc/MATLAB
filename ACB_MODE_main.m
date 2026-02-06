% ACB_MODE_main.m
% Adaptive-Constrained Bi-objective Multi-objective DE for VRP-C

clear; clc; rng(0);

%% ========================= PARAMETERS =========================
instance_name = 'rc102.txt';   % input instance (Solomon style)
NP = 50;                       % population size
GenMax = 100;                  % max generations
F0 = 1.0; Cr0 = 0.2;           % base DE parameters (mutate/crossover)
pLS = 0.5;                     % probability to run local search on offspring
archiveMax = 200;              % archive capacity
TR_default = 15;               % default TR for LNS (can be adapted)

% Vehicle/driver parameters (example for RC1)
v_num_full = [10, 40];         % [min, max] full-time drivers
v_num_part = [1, 15];          % [min, max] part-time drivers
cap_full = 200; cap_part = 150;
rate_full = 2; rate_part = 3;
r2_part = 100; e2_part = 200;   % part-time window (RC1 example)

% ACH (adaptive constraint balance) parameters (paper-style heuristics)
ACH.k = 5;                     % scale factor in adjustment (paper suggests small ints)
ACH.alpha0 = 0.9;              % initial epsilon multiplier
ACH.decay = 0.98;              % per-gen decay factor (if applicable)
ACH.eps = 1.0;                 % initial epsilon threshold (updated adaptively)
ACH.eps_min = 1e-6;

% other parameters
maxLSiter = 5;                 % max local search improvements per triggered solution
verbose = true;

%% ========================= 1. Read data (Solomon format) =========================
data = importdata(instance_name);
center = data(1,2:3);          % depot coords
E_center = data(1,5); L_center = data(1,6);
customers = data(2:end,2:3);
cusnum = size(customers,1);
demands = data(2:end,4);
a = data(2:end,5);
b = data(2:end,6);
s = data(2:end,7);             % service durations
vertexs = [center; customers];
dist = squareform(pdist(vertexs)); % Euclidean distance matrix

%% ========================= 2. Initialization (three-stage) =========================
% We'll generate initial population by the three-stage approach described in paper:
% (1) time-window / capacity grouping -> partition customers into full/part candidate sets
% (2) geo/time clustering inside groups (kmeans)
% (3) greedy/improved-CW (with timewindow checks) to produce initial routes, then encode to 3xN form

population = cell(NP,1);
objs = zeros(NP,2);
vio = zeros(NP,1);

for i=1:NP
    % Stage A: time-window + capacity grouping
    [full_cus, part_cus] = grouping_stage(customers, demands, a, b, cap_part, r2_part, e2_part);
    % Stage B: clustering inside groups to get subgroups (kmeans on geographic & time)
    [full_clusters] = cluster_stage(full_cus, customers, a, 3);  % K can be tuned
    [part_clusters] = cluster_stage(part_cus, customers, a, 2);
    % Stage C: greedy/improved CW per cluster to form initial routes
    enc = threeStage_build_routes(full_clusters, part_clusters, dist, demands, a, b, s, cap_full, cap_part, v_num_full, v_num_part);
    % Small random perturbations to increase diversity
    enc = perturb_encoding(enc);
    population{i} = enc;
end

% compute initial objectives & build archive (single struct with sols/objs/vio)
[f1_init, f2_init, vio_init] = calObjs(population, v_num_full, v_num_part, cap_full, cap_part, ...
    rate_full, rate_part, r2_part, e2_part, demands, a, b, s, E_center, L_center, dist);
archive = struct();
archive.sols = population;
archive.objs = [f1_init, f2_init];
archive.vio = vio_init;

%% ========================= 3. Main ACB-MODE loop =========================
for gen = 1:GenMax
    % adaptive DE parameters
    F = max(0.5, min(1.5, F0 + 0.05*randn));
    CR = max(0.05, min(0.9, Cr0 + 0.05*randn));
    
    % compute feasible ratio
    rf = sum(archive.vio==0)/max(1,length(archive.vio));
    % update ACH epsilon adaptively
    ACH = updateACH(ACH, archive, gen, GenMax);
    
    new_pop = cell(NP,1);
    new_objs = zeros(NP,2);
    new_vio = zeros(NP,1);
    
    for i=1:NP
        % Neighborhood selection using R(p|q) concept
        neigh = neighborhood_by_R(archive, i, NP, rf);
        [p1,p2,p3] = select_parents_from_neighborhood(neigh, length(archive.sols));
        x1 = archive.sols{p1}; x2 = archive.sols{p2}; x3 = archive.sols{p3};
        xi = archive.sols{min(i,length(archive.sols))};
        % Hybrid DE in continuous mapping -> decode -> repair
        u = hybridDE_and_decode(x1,x2,x3,xi,F,CR,cusnum,v_num_full,v_num_part);
        % possibly local search (neighborhood-based LNS or 2-opt)
        if rand < pLS
            u = local_search_triggered(u, demands, a, b, s, cap_full, cap_part, TR_default, dist, ACH);
        end
        % ensure driver constraints (split/merge)
        u = ensureDriverConstraints(u, v_num_full, v_num_part);
        % evaluate
        [f1,f2,viol] = calObj(u, v_num_full, v_num_part, cap_full, cap_part, ...
            rate_full, rate_part, r2_part, e2_part, demands, a, b, s, E_center, L_center, dist);
        new_pop{i} = u; new_objs(i,:) = [f1,f2]; new_vio(i) = viol;
    end
    
    % Combine and selection via ACH-based comparator (epsilons + CDP integration)
    all_pop = [archive.sols(:); new_pop(:)];
    all_objs = [archive.objs; new_objs];
    all_vio = [archive.vio(:); new_vio];
    % selection: prefer feasible; else epsilon-based compare; then nondominated & diversity
    [next_pop, next_objs, next_vio] = selection_ACB_MODE(all_pop, all_objs, all_vio, NP, ACH);
    
    % update archive as single struct
    archive.sols = next_pop;
    archive.objs = next_objs;
    archive.vio = next_vio;
    
    if verbose
        feas = sum(next_vio==0);
        fprintf('Gen %d: feasible=%d, minF1=%.2f, minF2=%.2f, eps=%.4g\n', gen, feas, min(next_objs(:,1)), min(next_objs(:,2)), ACH.eps);
    end
end

%% ========================= 4. Results & Plotting =========================
% get feasible front
feas_idx = find(archive.vio==0);
if isempty(feas_idx)
    sols_show = archive.sols;
    objs_show = archive.objs;
else
    sols_show = archive.sols(feas_idx);
    objs_show = archive.objs(feas_idx,:);
end

hv = calculateHV(objs_show);
E_val = calculateE(objs_show);
fprintf('\nFinal HV=%.6f, E=%.6f\n', hv, E_val);

% draw best (min F1) solution
[~, idx_best] = min(objs_show(:,1));
best_sol = sols_show{idx_best};
figure; drawRoutes(best_sol, vertexs);
title('Result routes (full-time: blue solid; part-time: cyan dashed)');

% Pareto front scatter
figure;
scatter(archive.objs(:,1), archive.objs(:,2), 36, 'r', 'filled'); hold on;
xlabel('F1 (Total distance)'); ylabel('F2 (Total pay)');
title('Pareto candidate set');
grid on;

