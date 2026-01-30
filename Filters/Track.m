classdef (Abstract) Track < handle
    %Generic EKF-type track with weighted association base class

    properties (Abstract)
        x
        P
    end

    methods (Abstract)
        obj = Extrapolate(obj);
        d2 = getChi2Dist(obj, x, R);
        S = getS(obj,R);
        obj = WeightedUpdate(obj, Measurements, R, Weights);
    end
end