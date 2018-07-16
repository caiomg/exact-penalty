function s = correct_step_to_bounds(x, s0, bl, bu, bl_active, bu_active)

    if isempty(bl)
        bl = -inf(size(x));
    end
    if isempty(bu)
        bu = inf(size(x));
    end
    if nargin < 5 || isempty(bl_active)
       bl_active = false(size(x)); 
    end
    if nargin < 6 || isempty(bu_active)
       bu_active = false(size(x)); 
    end
    slack_u = bu - x;
    slack_l = bl - x;
    
    s = min(slack_u, max(slack_l, s0));
    s(bl_active) = slack_l(bl_active);
    s(bu_active) = slack_u(bu_active);
    
end
    