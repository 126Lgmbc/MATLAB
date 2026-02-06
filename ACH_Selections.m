% ACH_Selection.m
% Contains ACH-based selection and neighborhood mechanisms

%% ---------- ACH update ----------
function ACH = updateACH(ACH, archive, gen, GenMax)
    % adaptive epsilon update: reduce slowly, but if no feasible solutions, raise to min violation
    if isempty(archive.vio)
        ACH.eps = ACH.eps * ACH.decay;
        ACH.eps = max(ACH.eps_min, ACH.eps);
        return;
    end
    feasible_ratio = sum(archive.vio==0)/length(archive.vio);
    if feasible_ratio == 0
        % find minimal violation among population
        min_vio = min(archive.vio);
        ACH.eps = max(ACH.eps, min_vio * ACH.alpha0);
    else
        % reduce epsilon gradually when feasible solutions exist
        ACH.eps = ACH.eps * ACH.decay;
        ACH.eps = max(ACH.eps_min, ACH.eps);
        % if feasible ratio large and improvement observed, relax slightly to allow exploration
        if feasible_ratio > 0.5
            ACH.eps = ACH.eps * (1 + 0.01*randn);
        end
    end
end
%% ---------- selection using ACH (epsilon constraint handling + nondom + diversity) ----------
function [next_pop, next_objs, next_vio] = selection_ACB_MODE(all_pop, all_objs, all_vio, NP, ACH)
    % first partition feasible & infeasible
    feasible_idx = find(all_vio==0);
    infeasible_idx = find(all_vio>0);
    selected = [];
    % among feasible: nondominated sorting & crowding (simple quick)
    if ~isempty(feasible_idx)
        feas_objs = all_objs(feasible_idx,:);
        fronts = nondominated_fronts(feas_objs);
        for f=1:length(fronts)
            idxs = feasible_idx(fronts{f});
            % if adding all fit
            if length(selected)+length(idxs) <= NP
                selected = [selected; idxs];
            else
                % crowding distance selection
                remain = NP - length(selected);
                cds = crowding_distance(all_objs(idxs,:));
                [~, ord] = sort(cds,'descend');
                selected = [selected; idxs(ord(1:remain))];
                break;
            end
        end
    end
    % if not full, select from infeasible using epsilon-comparator (paper's epsilon CHT)
    if length(selected) < NP && ~isempty(infeasible_idx)
        % compute comparator value: prefer smaller vio first if < eps, else use weighted rank of objectives
        % Here implement sample rule: sort by (vio <= ACH.eps) first, then by vio, then by sum(obj)
        small_vio = infeasible_idx(all_vio(infeasible_idx) <= ACH.eps);
        big_vio = setdiff(infeasible_idx, small_vio);
        candidates = [small_vio; big_vio];
        needed = NP - length(selected);
        if length(candidates) <= needed
            selected = [selected; candidates];
        else
            % sort small_vio by obj sum asc
            if ~isempty(small_vio)
                sums = sum(all_objs(small_vio,:),2);
                [~, ord] = sort(sums);
                take = small_vio(ord);
                selected = [selected; take(1:min(needed,length(take)))];
            end
            if length(selected) < NP
                % fill remaining by smallest violation
                remainNeeded = NP - length(selected);
                [~, ord2] = sort(all_vio(big_vio));
                selected = [selected; big_vio(ord2(1:min(remainNeeded,length(ord2))))];
            end
        end
    end
    % fallback if still not enough
    if length(selected) < NP
        allidx = 1:length(all_pop);
        remain = setdiff(allidx, selected);
        add = remain(1:min(NP-length(selected), length(remain)));
        selected = [selected; add];
    end
    selected = selected(1:NP);
    next_pop = all_pop(selected);
    next_objs = all_objs(selected,:);
    next_vio = all_vio(selected);
end
%% ---------- neighborhood_by_R (R(p|q) idea) ----------
function idxs = neighborhood_by_R(archive, curIdx, NP, rf)
    nA = length(archive.sols);
    if nA == 0, idxs = []; return; end
    curIdx = min(max(1, curIdx), nA);
    curV = archive.vio(curIdx);
    Rvals = inf(nA,1);
    for j=1:nA
        if j == curIdx, Rvals(j) = -inf; continue; end
        if archive.vio(j) < curV
            Rvals(j) = 0;
        elseif archive.vio(j) == curV
            rp = sum(archive.objs(curIdx,:)); rq = sum(archive.objs(j,:));
            if rq == 0, Rvals(j) = 1e6; else Rvals(j) = rp/rq; end
        else
            Rvals(j) = Inf;
        end
    end
    candidates = find(isfinite(Rvals));
    if isempty(candidates)
        idxs = setdiff(1:nA, curIdx);
    else
        [~, ord] = sort(Rvals(candidates), 'descend');
        Np = min(max(5, round(NP*0.5)), length(candidates));
        sel = candidates(ord(1:Np));
        idxs = sel(:)';
    end
end

%% ---------- robust parent selection ----------
function [p1,p2,p3] = select_parents_from_neighborhood(neighIdx, archiveSize)
    n = length(neighIdx);
    if n == 0
        % fallback to random distinct from archive
        if archiveSize < 3
            p1 = 1; p2 = 1; p3 = 1;
        else
            rp = randperm(archiveSize,3);
            p1 = rp(1); p2 = rp(2); p3 = rp(3);
        end
        return;
    elseif n == 1
        p1 = neighIdx(1); p2 = neighIdx(1); p3 = neighIdx(1);
        return;
    elseif n == 2
        rp = neighIdx(randperm(2));
        p1 = rp(1); p2 = rp(2); p3 = rp(1);
        return;
    else
        rp = neighIdx(randperm(n,3));
        p1 = rp(1); p2 = rp(2); p3 = rp(3);
        return;
    end
end