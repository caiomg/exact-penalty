function t = infinity_tr_radius_breakpoint(d, radius, s0)

    assert (radius > 0)
    if nargin < 3
        s0 = zeros(size(d));
    elseif ~isempty(find(s0 > radius, 1))
        warning('cmg:out_radius', 'Initial point out of radius');
    end

    increasing_coordinates = d > 0;
    decreasing_coordinates = d < 0;
    
    inc_max = min((radius - s0(increasing_coordinates))./ ...
                  d(increasing_coordinates));
    dec_max = min(-(radius - s0(decreasing_coordinates))./ ...
                  d(decreasing_coordinates));
    
    t = min([inc_max, dec_max]);

    assert(isempty(t) || (numel(t) == 1 && t >= 0));
    
end