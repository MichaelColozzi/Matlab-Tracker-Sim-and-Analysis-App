classdef JPDA < Tracker
    % basic Joint Probabilistic Data Association tracker class

    properties
        MostRecentId
        Tracks
        PD
        Beta
        A
        Q
        H
        InitialCovPad
        LoggingCell
        CurrLogCell
        currentFrame
    end

    methods
        function obj = JPDA(PD,Beta,A,Q,H,InitialCovPad)
            %starts with no tracks
            obj.MostRecentId = 1;
            obj.Tracks = {};
            obj.PD = PD;
            obj.Beta = Beta;
            obj.A = A;
            obj.Q = Q;
            obj.H = H;
            obj.InitialCovPad = InitialCovPad; %Added to state covariance of new tracks upon initialization
            obj.LoggingCell = {};
            obj.CurrLogCell = cell(1,8);
            obj.currentFrame = 0;
        end

        function obj = Extrapolate(obj)
            if isempty(obj.Tracks)
                return
            end

            for track = obj.Tracks{:}'
                track.Extrapolate();
            end
        end

        function obj = UpdateWithMeasurements(obj,Measurements,R)
            
            %If no tracks exist, create new ones from measurements.
            if isempty(obj.Tracks)
                obj.Tracks = cell(size(Measurements));
                for n = 1:numel(Measurements)
                    pinvH = pinv(obj.H);
                    obj.Tracks{n} = IdealKFTrack(obj.MostRecentId,...
                                                 pinvH * Measurements{n},...
                                                 pinvH * R{n} * pinvH' + obj.InitialCovPad,...
                                                 obj.A,obj.Q,obj.H);
                    obj.MostRecentId = obj.MostRecentId + 1;
                end
                obj.LogMeas(Measurements,R,[]);
                obj.LogCurr();
                return
            end

            AssociationCostMatrix = zeros(numel(obj.Tracks),numel(Measurements));
            LogDetS = zeros(numel(obj.Tracks),numel(Measurements));
            
            i = 0;
            for track = obj.Tracks(:)'
                i = i + 1;
                j = 0;
                for meas = Measurements(:)'
                    j = j+1;
                    LogDetS(i,j) = log(det(track{1}.getS(R{j})));
                    AssociationCostMatrix(i,j) = track{1}.getChi2Dist(meas{1},R{j}) + LogDetS(i,j);
                end
            end
            mlgate = 2* log(obj.PD./((2*pi).^(size(obj.H,1)/2).*obj.Beta));
            AssociationCostMatrix((AssociationCostMatrix>mlgate)) = inf;
            [M,uR,uC] = matchpairs(AssociationCostMatrix,1e6,'min');
            W = exp(-0.5*sum(AssociationCostMatrix(sub2ind(...
                                                   size(AssociationCostMatrix),...
                                                   M(:,1),M(:,2))),"all")).*ones(size(M,1),1);
            MAug = [M,W];
            for assocPair = M'
                miniCost = AssociationCostMatrix;
                miniCost(assocPair') = inf;
                [Mmini,~,~] = matchpairs(miniCost,1e6,'min');
                Wmini = exp(-0.5*sum(AssociationCostMatrix(sub2ind(...
                                                   size(AssociationCostMatrix),...
                                                   Mmini(:,1),Mmini(:,2))),"all")).*ones(size(Mmini,1),1);
                MAug = [MAug;[Mmini,Wmini]];
            end
            for k=1:size(AssociationCostMatrix,1)
                MAug(MAug(:,1)==k,3)=MAug(MAug(:,1)==k,3)./sum(MAug(MAug(:,1)==k,3));
            end
            assocweights = zeros(size(AssociationCostMatrix));
            %for weightedpair = MAug'
            %    assocweights(weightedpair(1),weightedpair(2)) = assocweights(weightedpair(1),weightedpair(2)) + weightedpair(3);
            %end
            MAug(isnan(MAug(:,3)),:)=[];
            %Logging
            obj.LogMeas(Measurements,R,MAug(:,1:2));

            %unmatched measurements
            obj.LogBirths([length(obj.Tracks)+1:length(obj.Tracks)+length(uC)]);
            for unmatchedMeas = uC(:)'
                pinvH = pinv(obj.H);
                obj.Tracks{end+1} = IdealKFTrack(obj.MostRecentId, ...
                                                 pinvH * Measurements{unmatchedMeas},...
                                                 pinvH * R{unmatchedMeas} * pinvH' + obj.InitialCovPad,...
                                                 obj.A,obj.Q,obj.H);
                obj.MostRecentId = obj.MostRecentId + 1;
            end

            %unmatched tracks
            for umatchedTrk = uR(:)'
                obj.Tracks{umatchedTrk}.Extrapolate();
            end

            %Associations
            if ~isempty(AssociationCostMatrix)
            for k=1:size(AssociationCostMatrix,1)
                obj.Tracks{k}.WeightedUpdate(Measurements(MAug(MAug(:,1)==k,2)),...
                                                        R(MAug(MAug(:,1)==k,2)),...
                                                        MAug(MAug(:,1)==k,3));
            end
            end

            tidx = 0;
            tracks2delete = [];
            for track = obj.Tracks(:)'
                tidx = tidx + 1;
                if ((1-obj.PD)^(track{1}.PropogationsSinceLastUpdate))<(0.01/(1+track{1}.TrackLifetime))
                    tracks2delete(end+1) = tidx;
                end
            end
            obj.LogDeaths(tracks2delete);
            obj.Tracks(tracks2delete) = [];

            obj.LogCurr();
        end
        function obj = LogMeas(obj,Measurements,R,Associations)
            %Trk id's, trk states, trk covariances, measurements,
            %measurement covariances, Associations, births, deaths
            Ids = getTrkIds(obj);
            [xs,Ps] = obj.getTrkStates();
            meas =cat(2,Measurements{:});
            measCov = cat(3,R{:});
            obj.CurrLogCell(1:6) = {Ids,xs,Ps,meas,measCov,Associations};
        end
        function obj = LogBirths(obj,BirthIndices)
            %Trk id's, trk states, trk covariances, measurements,
            %measurement covariances, Associations, births, deaths
            obj.CurrLogCell(7) = {BirthIndices};
        end
        function obj = LogDeaths(obj,DeathIndices)
            %Trk id's, trk states, trk covariances, measurements,
            %measurement covariances, Associations, births, deaths
            obj.CurrLogCell(8) = {DeathIndices};
        end
        function obj = LogCurr(obj)
            if ~isempty(obj.LoggingCell)
                obj.LoggingCell(end+1,:) = obj.CurrLogCell;
            else
                obj.LoggingCell = obj.CurrLogCell;
            end
        end
        function Ids = getTrkIds(obj)
            Ids = zeros(size(obj.Tracks));
            for n = 1:length(obj.Tracks)
                Ids(n) = obj.Tracks{n}.id;
            end
        end
        function [xs,Ps] = getTrkStates(obj)
            xs = zeros(4,length(obj.Tracks));
            Ps = zeros(4,4,length(obj.Tracks));
            for n = 1:length(obj.Tracks)
                xs(:,n) = obj.Tracks{n}.x;
                Ps(:,:,n) = obj.Tracks{n}.P;
            end
        end
    end
end