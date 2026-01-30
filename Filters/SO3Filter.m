classdef SO3Filter < Track
    %attitude filter

    properties
        DCM
        logP
        Omega
        logDCM
        AOmega
        dt
        Q
        H
    end

    methods
        function obj = SO3Filter(DCM,...
                                 logP,...
                                 Omega,...
                                 AOmega,...
                                 dt, ...
                                 Q,...
                                 H)
            %provide values for all the major filter parameters
            obj.DCM = DCM;
            obj.logP = logP;
            obj.Omega = Omega;
            obj.logDCM = logm(DCM);
            obj.AOmega = AOmega;
            obj.dt = dt;
            obj.Q = Q;
            obj.H = H;
        end

        function obj = WeightedUpdate(obj, Measurements, R, Weights)

            obj.Extrapolate();
            if isempty(Measurements);return ;end

            Beta0 = 1-sum(Weights);
            innovations = cell(size(Measurements));
            KalmanGains = cell(size(Measurements));
            InnovationCovariances = cell(size(Measurements));
            SingleMeasx = cell(size(Measurements));
            SingleMeasP = cell(size(Measurements));
            
            UpdatedState = obj.x;
            for n = 1:numel(Measurements) 
                innovations{n} = Measurements{n} - obj.H * obj.x;
                InnovationCovariances{n} = obj.H * obj.P * obj.H' + R{n};
                KalmanGains{n} = obj.P * obj.H' * (InnovationCovariances{n})^-1;
                correction = KalmanGains{n} * innovations{n};
                SingleMeasx{n} = obj.x + correction;
                UpdatedState = UpdatedState + Weights(n) .* correction;
                temp1 = KalmanGains{n}*obj.H;
                temp2 = eye(size(temp1,1)) - temp1;
                SingleMeasP{n} = temp2*obj.P*temp2' + KalmanGains{n} * R{n} * KalmanGains{n}';
            end
            totalCorrection = obj.x - UpdatedState;
            Pnew = Beta0 .* (obj.P + totalCorrection * totalCorrection');
            for n = 1:numel(Measurements) 
                correction = SingleMeasx{n}-UpdatedState;
                Pnew = Pnew + Weights(n) .* (SingleMeasP{n}+correction*correction');
            end
            obj.x = UpdatedState;
            obj.P = Pnew;
        end
        function obj = Extrapolate(obj)
            obj.DCM = expm(obj.dt*obj.Omega)*obj.DCM;
            logDCMNew = logm(obj.DCM);
            Aeff = skew(logDCMNew)*pinv(skew(obj.logDCM));
            obj.Omega = skew(obj.AOmega * skew(obj.Omega));
        end

        function w = skew(O)
            if all(size(O)==[3,3])
                w=[-O(2,3);O(1,3);-O(1,2)];
            else
                w = zeros(3);
                w(1,2) = -O(3);
                w(1,3) = O(2);
                w(2,3) = -O(1);
                w = w - w';
            end
        end

    end
end