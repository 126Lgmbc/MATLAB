% SolutionOperators.m
% Contains all solution generation and manipulation functions

%% ---------- Grouping stage ----------
function [full_cus, part_cus] = grouping_stage(customers, demands, a, b, cap_part, r2_part, e2_part)
    % select clients eligible for part-time by time-window inclusion and capacity heuristic
    N = size(customers,1);
    idx_part = find(a >= r2_part & b <= e2_part);  % strict inclusion
    if isempty(idx_part)
        % relax to overlap criterion
        idx_part = find(~(b < r2_part | a > e2_part));
    end
    % capacity heuristic: pick subset for part-time that cumulative demand per route plausible
    % For init, just keep idx_part but limit size to ensure min customers for full drivers
    part_cus = idx_part(:)';
    full_cus = setdiff(1:N, part_cus);
    % ensure minimal split
    if isempty(full_cus)
        % assign some to full randomly
        k = max(1, floor(0.2*N));
        full_cus = part_cus(1:k);
        part_cus(1:k) = [];
    end
end

%% ---------- Clustering stage ----------
function clusters = cluster_stage(cus_list, customers, a, K)
    % do a small kmeans on coords + time-start as feature (if too few customers, return single cluster)
    if isempty(cus_list)
        clusters = {};
        return;
    end
    if length(cus_list) <= K
        clusters = cell(1,length(cus_list));
        for i=1:length(cus_list), clusters{i} = cus_list(i); end
        return;
    end
    feats = [customers(cus_list, :), a(cus_list)];
    try
        idxk = kmeans(feats, K, 'MaxIter',100,'Replicates',3);
    catch
        idxk = randi(K, size(feats,1),1);
    end
    clusters = cell(1,K);
    for k=1:K
        clusters{k} = cus_list(idxk==k);
        if isempty(clusters{k}), clusters{k} = []; end
    end
    clusters = clusters(~cellfun(@isempty,clusters));
end

%% ---------- Build routes (3-stage combined) ----------
function enc = threeStage_build_routes(full_clusters, part_clusters, dist, demands, a, b, s, cap_full, cap_part, v_num_full, v_num_part)
    % For each cluster, run improved CW with time-window check, then assemble into global encoding
    all_cus = unique([cell2mat(full_clusters(:)'), cell2mat(part_clusters(:)')]);
    max_idx = max([all_cus, 0]);
    if max_idx == 0
        enc = zeros(3,0); return;
    end
    type = zeros(1,max_idx); v_id = zeros(1,max_idx); seq = zeros(1,max_idx);
    next_vid = 1;
    % full clusters -> closed routes
    for k=1:length(full_clusters)
        members = full_clusters{k}(:)';
        if isempty(members), continue; end
        routes = improvedCW_cluster(members, dist, demands, a, b, s, cap_full, 0);
        for r=1:length(routes)
            mem = routes{r};
            v_id(mem) = next_vid;
            type(mem) = 0;
            [~, ord] = sort(a(mem));
            seq(mem) = 1:length(mem);
            next_vid = next_vid + 1;
        end
    end
    % part clusters -> open routes
    for k=1:length(part_clusters)
        members = part_clusters{k}(:)';
        if isempty(members), continue; end
        routes = improvedCW_cluster(members, dist, demands, a, b, s, cap_part, 1);
        for r=1:length(routes)
            mem = routes{r};
            v_id(mem) = next_vid;
            type(mem) = 1;
            [~, ord] = sort(a(mem));
            seq(mem) = 1:length(mem);
            next_vid = next_vid + 1;
        end
    end
    enc = zeros(3, max_idx);
    enc(1,:) = type; enc(2,:) = v_id; enc(3,:) = seq;
end

function routes = improvedCW_cluster(customers, dist, demands, a, b, s, cap, routeOpen)
    % simplified CW saving per cluster with time-window feasibility check
    customers = customers(:)';
    n = length(customers);
    if n==0, routes={}; return; end
    savings = [];
    for i=1:n
        for j=i+1:n
            ci = customers(i); cj = customers(j);
            sij = dist(1,ci+1) + dist(1,cj+1) - dist(ci+1,cj+1);
            savings = [savings; i, j, sij];
        end
    end
    savings = sortrows(savings, -3);
    % init routes: each customer its own route
    routes = cell(1,n);
    for i=1:n, routes{i} = customers(i); end
    % merge greedily while capacity allows and rough time-window check
    idx = 1;
    while idx <= size(savings,1)
        i = savings(idx,1); j = savings(idx,2);
        ci = customers(i); cj = customers(j);
        ri = find(cellfun(@(r) any(r==ci), routes), 1);
        rj = find(cellfun(@(r) any(r==cj), routes), 1);
        if isempty(ri) || isempty(rj) || ri==rj
            idx = idx + 1; continue;
        end
        merged = [routes{ri}, routes{rj}];
        if sum(demands(merged)) <= cap
            % sort by earliest start time to try to satisfy time windows
            [~, ord] = sort(a(merged));
            cand = merged(ord);
            if checkRouteTimeFeasible(cand, dist, a, b, s, routeOpen)
                routes{ri} = cand;
                routes(rj) = [];
            end
        end
        idx = idx + 1;
    end
    routes = routes(~cellfun(@isempty, routes));
end

function enc = perturb_encoding(enc)
    % Robust perturbation: small inter-route swaps or intra-route random reorder
    if isempty(enc), return; end
    type = enc(1,:); v_id = enc(2,:); seq = enc(3,:);
    all_vid = unique(v_id); all_vid(all_vid==0) = [];
    if length(all_vid) >= 2
        if rand < 0.5
            % swap two customers belonging to two different vehicles
            vsel = randsample(all_vid,2);
            c1 = find(v_id==vsel(1)); c2 = find(v_id==vsel(2));
            if ~isempty(c1) && ~isempty(c2)
                idx1 = c1(randi(length(c1))); idx2 = c2(randi(length(c2)));
                tmpv = v_id(idx1); v_id(idx1) = v_id(idx2); v_id(idx2) = tmpv;
                % optionally swap seq as well
                tmps = seq(idx1); seq(idx1) = seq(idx2); seq(idx2) = tmps;
            end
        else
            % randomly shuffle sequence inside a randomly chosen vehicle
            v = all_vid(randi(length(all_vid)));
            cs = find(v_id==v);
            if length(cs) >= 2
                % ensure cs is a column vector index list
                cs = cs(:)';
                % take current sequence values
                tmp = seq(cs);
                % random permutation of indices 1..numel(cs)
                perm = randperm(numel(cs));
                % assign back permuted values
                seq(cs) = tmp(perm);
            end
        end
    end
    enc(2,:) = v_id; enc(3,:) = seq;
end

%% ---------- Hybrid DE on continuous mapped encoding + decode ----------
function u = hybridDE_and_decode(x1,x2,x3,xi,F,CR,cusnum,vnumFull,vnumPart)
    % Map 3xN discrete encoding -> continuous vector of same size in [0,1]
    % For simplicity we create cont matrix of same size: row1:type weight in [0,1], row2:vid weight, row3:seq weight
    D = size(x1,2);
    cont1 = map_encoding_to_cont(x1, D);
    cont2 = map_encoding_to_cont(x2, D);
    cont3 = map_encoding_to_cont(x3, D);
    % DE mutation
    v = cont1 + F*(cont2 - cont3);
    % crossover
    u_cont = map_encoding_to_cont(xi, D);
    jrand = randi(D);
    for j=1:D
        if rand < CR || j==jrand
            u_cont(:,j) = v(:,j);
        end
    end
    % decode continuous to discrete encoding
    u = decode_cont_to_encoding(u_cont, cusnum, vnumFull, vnumPart);
end

function cont = map_encoding_to_cont(enc, D)
    % enc is 3xN (may be shorter than D)
    % produce 3xD continuous ~[0,1]
    cont = zeros(3,D);
    if isempty(enc)
        cont = rand(3,D);
        return;
    end
    % row1: type in {0,1} -> map 0->0.1, 1->0.9 with small jitter
    cont(1,1:size(enc,2)) = enc(1,1:size(enc,2))*0.8 + 0.1 + 0.01*randn(1,size(enc,2));
    % row2: v_id -> normalized by max vid
    maxv = max(1,max(enc(2,:)));
    cont(2,1:size(enc,2)) = enc(2,1:size(enc,2))/maxv;
    % row3: seq -> normalized by max seq
    maxs = max(1,max(enc(3,:)));
    cont(3,1:size(enc,2)) = enc(3,1:size(enc,2))/maxs;
    % fill rest random if D > length(enc)
    if D > size(enc,2)
        cont(:, size(enc,2)+1:D) = rand(3, D - size(enc,2));
    end
    % clip
    cont = min(1, max(0, cont));
end

function enc = decode_cont_to_encoding(u_cont, cusnum, vnumFull, vnumPart)
    % decode cont into discrete 3xN encoding for N=cusnum
    D = size(u_cont,2);
    N = cusnum;
    enc = zeros(3,N);
    % if D < N: pad with random
    if D < N
        u_cont = [u_cont, rand(3, N - D)];
    end
    % row1: type threshold 0.5
    types = double(u_cont(1,1:N) > 0.5);
    % row2: v_id from cont2: we map to integer in [1, maxVidEstimate]
    maxVidEstimate = max(5, round(0.5 * N));
    vids = max(1, round(u_cont(2,1:N) * maxVidEstimate));
    % row3: seq -> ranks of cont value to produce order within vehicle; but initial seq assign as rank across customers; will be normalized later
    seq_raw = u_cont(3,1:N);
    % create initial enc
    enc(1,:) = types;
    enc(2,:) = vids;
    % produce seq per vehicle: sort by seq_raw inside each vid
    uniqv = unique(vids);
    for v = uniqv
        idxs = find(vids == v);
        [~, ord] = sort(seq_raw(idxs));
        seqvals = 1:length(idxs);
        enc(3, idxs(ord)) = seqvals;
    end
    % ensure seq nonzero
    enc(3, enc(3,:) == 0) = 1;
    % convert types to 0/1
    enc(1,:) = types;
end

%% ---------- Local search triggered (LNS + neighbor operators) ----------
function sol = local_search_triggered(sol, demands, a, b, s, cap_full, cap_part, TR, dist, ACH)
    % decide whether to run LNS or smaller ops: in paper it's triggered under conditions
    % Here we always run LNS with moderate TR but limited iterations
    if isempty(sol), return; end
    type = sol(1,:); v_id = sol(2,:); seq = sol(3,:);
    N = length(type);
    % compute correlation-based relevance matrix (simple): use distance + time overlap
    rel = compute_relevance(N, dist, a, b);
    % select TR nodes with highest aggregated relevance (sum over rows)
    agg = sum(rel,2);
    [~, ridx] = sort(agg,'descend');
    to_remove = ridx(1:min(TR, N));
    % remove them (unassign)
    v_id(to_remove) = 0; seq(to_remove) = 0;
    unassigned = to_remove;
    % repair by best insertion (consider type compatibility and time feasibility)
    % attempt to insert into existing vehicles (same type first), else create new
    all_vids = unique(sol(2,:)); all_vids(all_vids==0) = [];
    attempts = 0;
    while ~isempty(unassigned) && attempts < 5*length(unassigned)
        attempts = attempts + 1;
        % greedy choose the best (min cost increase) insertion among all unassigned
        best_gain = inf; best_c = -1; best_vid = -1; best_pos = -1; best_type = -1;
        for c = unassigned
            % try same-type vehicles first
            cand_types = [0,1];
            for t = cand_types
                vids_t = unique(sol(2, sol(1,:)==t)); vids_t(vids_t==0) = [];
                if isempty(vids_t)
                    % consider create new vehicle of type t
                    vid_try = max([sol(2,:), 0]) + 1;
                    pos_try = 1;
                    inc = estimate_insertion_cost_new(vid_try, c, sol, dist);
                    if inc < best_gain
                        best_gain = inc; best_c = c; best_vid = vid_try; best_pos = pos_try; best_type = t;
                    end
                else
                    for vv = vids_t
                        cs = find(sol(2,:)==vv);
                        for pos = 0:length(cs)
                            inc = estimate_insertion_cost(sol, vv, pos, c, dist);
                            if inc < best_gain
                                best_gain = inc; best_c = c; best_vid = vv; best_pos = pos; best_type = t;
                            end
                        end
                    end
                end
            end
        end
        if best_c == -1
            break;
        end
        % perform insertion
        sol(1,best_c) = best_type;
        sol(2,best_c) = best_vid;
        sol(3,best_c) = 9999; % will normalize later
        % remove from unassigned
        unassigned(unassigned==best_c) = [];
    end
    % normalize seq inside vehicles
    vids = unique(sol(2,:)); vids(vids==0) = [];
    for vv = vids
        cs = find(sol(2,:)==vv);
        if isempty(cs), continue; end
        [~, ord] = sort(a(cs));
        sol(3, cs(ord)) = 1:length(cs);
    end
end

function rel = compute_relevance(N, dist, a, b)
    % simple relevance: inverse distance + time window overlap score
    rel = zeros(N);
    for i=1:N
        for j=1:N
            if i==j, rel(i,j)=0; continue; end
            d = dist(i+1, j+1);
            overlap = max(0, min(b(i), b(j)) - max(a(i), a(j)));
            rel(i,j) = 1/(1+d) + overlap;
        end
    end
end

function inc = estimate_insertion_cost(sol, vid, pos, c, dist)
    % estimate distance increase if insert c into vehicle vid at position pos
    cs = find(sol(2,:)==vid);
    if isempty(cs)
        inc = dist(1, c+1); return;
    end
    if pos==0
        after = cs(1);
        inc = dist(1, c+1) + dist(c+1, after+1) - dist(1, after+1);
    elseif pos==length(cs)
        before = cs(end);
        inc = dist(before+1, c+1);
    else
        before = cs(pos); after = cs(pos+1);
        inc = dist(before+1, c+1) + dist(c+1, after+1) - dist(before+1, after+1);
    end
end

function inc = estimate_insertion_cost_new(vid, c, sol, dist)
    % cost of creating new vehicle and assigning c
    inc = dist(1, c+1);
end