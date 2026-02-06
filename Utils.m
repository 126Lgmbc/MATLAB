%% ---------- Non-dominated front for a set of 2D objs ----------
function fronts = nondominated_fronts(objs)
    N = size(objs,1);
    remaining = 1:N;
    fronts = {};
    while ~isempty(remaining)
        front = [];
        for i = remaining
            dominated = false;
            for j = remaining
                if i==j, continue; end
                if (objs(j,1) <= objs(i,1) && objs(j,2) <= objs(i,2)) && (objs(j,1) < objs(i,1) || objs(j,2) < objs(i,2))
                    dominated = true; break;
                end
            end
            if ~dominated, front = [front, i]; end
        end
        fronts{end+1} = front;
        remaining = setdiff(remaining, front);
    end
end

function cds = crowding_distance(objs)
    % simple crowding distance for 2D objs (higher means more isolated)
    N = size(objs,1);
    cds = zeros(N,1);
    for m = 1:2
        [vals, idx] = sort(objs(:,m));
        cds(idx(1)) = cds(idx(1)) + Inf;
        cds(idx(end)) = cds(idx(end)) + Inf;
        rngs = vals(end) - vals(1);
        if rngs == 0, rngs = 1; end
        for i = 2:N-1
            cds(idx(i)) = cds(idx(i)) + (vals(i+1) - vals(i-1))/rngs;
        end
    end
end
%% ---------- drawRoutes (visualization) ----------
function drawRoutes(sol, vertexs)
    type = sol(1,:); v_id = sol(2,:); seq = sol(3,:);
    center = vertexs(1,:); customers = vertexs(2:end,:);
    figure; hold on;
    plot(center(1), center(2), 'ks', 'MarkerSize',10, 'MarkerFaceColor',[0.3 0.3 0.3]);
    full_cus = find(type==0);
    plot(customers(full_cus,1), customers(full_cus,2), 'ko', 'MarkerFaceColor',[0.5 0.5 0.5]);
    part_cus = find(type==1);
    plot(customers(part_cus,1), customers(part_cus,2), 'ks', 'MarkerFaceColor',[0.8 0.8 0.8]);
    full_v = unique(v_id(type==0)); full_v(full_v==0) = [];
    % #0012c3
    full_color = [0, 18, 195] / 255;
    for ii=1:length(full_v)
        vid = full_v(ii);
        cus = find(type==0 & v_id==vid);
        [~, idx] = sort(seq(cus)); route = cus(idx);
        coords = [center; customers(route,:); center];
        plot(coords(:,1), coords(:,2), 'Color', full_color, 'LineWidth',1.5);
    end
    part_v = unique(v_id(type==1)); part_v(part_v==0) = [];
    % #ff84e8
    part_color = [255, 132, 232] / 255;
    for ii=1:length(part_v)
        vid = part_v(ii);
        cus = find(type==1 & v_id==vid);
        [~, idx] = sort(seq(cus)); route = cus(idx);
        coords = [center; customers(route,:)];
        plot(coords(:,1), coords(:,2), 'Color', part_color, 'LineWidth',1.5, 'LineStyle','--');
    end
    legend('depot','full customers','part customers','full routes','part routes');
    hold off; grid on;
end