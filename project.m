function zs = project(ps, ys) 
    % Finds the closest point on the boundary to the points in ys.
    
    nps = size(ps, 2);
    nys = size(ys, 2);
    
    as = ps(:,1:nps-1);
    bs = ps(:,2:nps) - as;
    
    zs = zeros(2, nys);
    
    for i = 1:nys
        y = ys(:,i);
        yas = y - as;
        ts = sum(yas .* bs, 1) ./ sum(bs .* bs, 1);
        ts = max(min(ts, 1), 0);
        
        yzs = as + ts .* bs;
        
        dists = vecnorm(yzs - y, 2);
        [~,index] = min(dists);
        zs(:, i) = yzs(:,index);
    end            
end