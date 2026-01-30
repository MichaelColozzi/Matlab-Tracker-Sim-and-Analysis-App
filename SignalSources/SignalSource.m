classdef (Abstract) SignalSource < handle

    methods (Abstract)
        obj = Update(obj);
        signalOuput = getSignalOutput(obj);
    end
end