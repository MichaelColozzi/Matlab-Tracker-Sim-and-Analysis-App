classdef Markov < handle
    %wrapper class which implements a markov chain of different signal
    %modes

    properties
        CurrentState
        transitionMatrix
        SignalModes
    end

    methods
        function obj = Markov(transitionMat,initialState,SignalModes)
            %Markov Construct an instance of this class
            %   Detailed explanation goes here
            obj.CurrentState = initialState;
            obj.transitionMatrix = transitionMat;
            obj.SignalModes = SignalModes;
        end

        function obj = Update(obj)
            obj.CurrentState = randsample(1:length(obj.transitionMatrix),1,true,obj.transitionMatrix(obj.CurrentState,:));
            obj.SignalModes{obj.CurrentState}.Update();
        end
        function signal = getSignalOutput(obj)
            signal = obj.SignalModes{obj.CurrentState}.getSignalOutput();
        end
    end
end