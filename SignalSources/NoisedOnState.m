classdef NoisedOnState < SignalSource
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Distribution
        value
    end

    methods
        function obj = NoisedOnState(Distribution)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Distribution = Distribution;
        end

        function obj = Update(obj)
            obj.value = random(obj.Distribution);
        end
        function output = getSignalOutput(obj)
            output = obj.value;
        end
    end
end