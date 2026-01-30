classdef (Abstract) Tracker < handle
    %Tracker Classes inherit from this

    properties (Abstract)
        Tracks
    end

    methods (Abstract)
        obj = Extrapolate(obj);
        obj = UpdateWithMeasurements(obj, Measurements, R);
    end
end