classdef OffState < SignalSource
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods
        function obj = OffState()
        end

        function obj = Update(obj)
            obj = obj;
        end
        function output = getSignalOutput(obj)
            output = 0;
        end
    end
end