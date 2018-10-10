classdef tr_model < handle
    %TR_MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        points_abs
        fvalues
        points_shifted
        cached_points
        cached_fvalues
        tr_center
        radius
        pivot_polynomials
        pivot_absvalues
        modeling_polynomials
    end
    
    methods
        function self = tr_model(points, fvalues, radius)
            self.points_abs = points;
            self.fvalues = fvalues;
            self.radius = radius;
            self.tr_center = 1;
        end
    end
    
end

