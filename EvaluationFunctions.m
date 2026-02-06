% EvaluationFunctions.m
% Contains objective evaluation and constraint checking functions

%% ---------- calObj (single solution) ----------
function [f1, f2, vio] = calObj(sol, v_num_full, v_num_part, cap_full, cap_part, ...
    rate_full, rate_part, r2_part, e2_part, demands, a, b, s, E, L, dist)
    if isempty(sol)
        f1 = Inf; f2 = Inf; vio = Inf; return;
    end
    type = sol(1,:); v_id = sol(2,:); seq = sol(3,:);
    f1 = 0; f2 = 0; vio = 0;
    full_v = unique(v_id(type==0)); full_v(full_v==0) = [];
    part_v = unique(v_id(type==1)); part_v(part_v==0) = [];
    vio = vio + max(0, v_num_full(1) - length(full_v)) * 100;
    vio = vio + max(0, length(full_v) - v_num_full(2)) * 100;
    vio = vio + max(0, v_num_part(1) - length(part_v)) * 100;
    vio = vio + max(0, length(part_v) - v_num_part(2)) * 100;
    % full routes
    for vv = full_v
        cus = find(type==0 & v_id==vv);
        if isempty(cus), continue; end
        [~, idx] = sort(seq(cus));
        route = cus(idx);
        % distance (closed)
        route_dist = dist(1, route(1)+1);
        for k=2:length(route)
            route_dist = route_dist + dist(route(k-1)+1, route(k)+1);
        end
        route_dist = route_dist + dist(route(end)+1, 1);
        f1 = f1 + route_dist;
        [time_v, cap_v, work_time] = checkConstraints_route(route, 0, cap_full, demands, a, b, s, E, L, r2_part, e2_part, dist);
        vio = vio + time_v + cap_v;
        f2 = f2 + rate_full * work_time;
    end
    % part routes
    for vv = part_v
        cus = find(type==1 & v_id==vv);
        if isempty(cus), continue; end
        [~, idx] = sort(seq(cus));
        route = cus(idx);
        route_dist = dist(1, route(1)+1);
        for k=2:length(route)
            route_dist = route_dist + dist(route(k-1)+1, route(k)+1);
        end
        f1 = f1 + route_dist;
        [time_v, cap_v, work_time] = checkConstraints_route(route, 1, cap_part, demands, a, b, s, E, L, r2_part, e2_part, dist);
        vio = vio + time_v + cap_v;
        f2 = f2 + rate_part * work_time;
    end
end

function [time_vio, cap_vio, work_time] = checkConstraints_route(route, type, cap, demands, a, b, s, E, L, r2, e2, dist)
    time_vio = 0;
    cap_vio = max(0, sum(demands(route)) - cap);
    if isempty(route), work_time = 0; return; end
    arr_time = E + dist(1, route(1)+1);
    if type == 1
        arr_time = max(arr_time, r2);
        time_vio = time_vio + max(0, arr_time - e2);
    end
    service_end = max(arr_time, a(route(1))) + s(route(1));
    time_vio = time_vio + max(0, a(route(1)) - arr_time) + max(0, service_end - b(route(1)));
    for i=2:length(route)
        arr_time = service_end + dist(route(i-1)+1, route(i)+1);
        if type == 1
            time_vio = time_vio + max(0, arr_time - e2);
        end
        service_end = max(arr_time, a(route(i))) + s(route(i));
        time_vio = time_vio + max(0, a(route(i)) - arr_time) + max(0, service_end - b(route(i)));
    end
    if type == 0
        return_time = service_end + dist(route(end)+1, 1);
        time_vio = time_vio + max(0, return_time - L);
        work_time = return_time - E;
    else
        work_time = service_end - max(E, r2);
    end
end

%% ---------- calObjs batch ----------
function [f1, f2, vio] = calObjs(pop, v_num_full, v_num_part, cap_full, cap_part, ...
    rate_full, rate_part, r2_part, e2_part, demands, a, b, s, E, L, dist)
    NP = numel(pop);
    f1 = zeros(NP,1); f2 = zeros(NP,1); vio = zeros(NP,1);
    for i=1:NP
        [f1(i), f2(i), vio(i)] = calObj(pop{i}, v_num_full, v_num_part, cap_full, cap_part, ...
            rate_full, rate_part, r2_part, e2_part, demands, a, b, s, E, L, dist);
    end
end

%% ---------- HV and E indicators ----------
function hv = calculateHV(objs)
    if isempty(objs), hv = 0; return; end
    ref = max(objs) + 1;
    hv = 0;
    for i=1:size(objs,1)
        hv = hv + (ref(1)-objs(i,1))*(ref(2)-objs(i,2));
    end
    hv = hv / (ref(1)*ref(2));
end

function e = calculateE(objs)
    if size(objs,1)<2, e=0; return; end
    e = 0;
    for m=1:2
        f = objs(:,m);
        f_norm = (f - min(f))/(max(f)-min(f)+eps);
        e = e + std(f_norm);
    end
end

function ok = checkRouteTimeFeasible(route, dist, a, b, s, routeTypeOpen)
% CHECKROUTETIMEFEASIBLE Checks the temporal feasibility of a single route.
%
% This function performs a rapid feasibility check for a given sequence of 
% customers, assuming the departure time from the depot is 0.
%
% Inputs:
%   route         - Vector of customer indices (1..N) in sequence.
%   dist          - Distance matrix (index 1 is depot; customer i is index i+1).
%   a, b, s       - Vectors for time window start, end, and service duration.
%   routeTypeOpen - 0 for closed routes (Full-time); 1 for open routes (Part-time).
%
% Output:
%   ok - Boolean: true if the sequence is temporally feasible; false otherwise.
%
% Details:
%   - Uses a greedy approach to calculate the earliest arrival time at each node.
%   - Checks if the arrival time exceeds the customer's upper bound (b).
%   - Accommodates service time (s) and waiting time (max(arrival, a)).

    ok = true;
    if isempty(route)
        return;
    end

    % Ensure route is a row vector
    route = route(:)';
    
    % Initialize time (assuming departure from depot at time 0)
    t = 0;
    
    % Travel to the first customer
    firstCust = route(1);
    t = t + dist(1, firstCust + 1);
    
    % Check if arrival time is within the time window's upper bound
    if t > b(firstCust)
        ok = false; 
        return;
    end
    
    % Start service at max(arrival, ready_time) and add service duration
    t = max(t, a(firstCust)) + s(firstCust);
    
    % Process subsequent customers in the sequence
    for k = 2:length(route)
        prev = route(k-1);
        cur  = route(k);
        t = t + dist(prev + 1, cur + 1);
        
        % Feasibility check
        if t > b(cur)
            ok = false; 
            return;
        end
        
        % Update current time after service
        t = max(t, a(cur)) + s(cur);
    end
    
    % Optional: Check return to depot for closed routes (Full-time)
    if routeTypeOpen == 0
        last = route(end);
        t = t + dist(last + 1, 1);
        % Note: If a global time limit (L) is defined, it should be checked here.
    end
end