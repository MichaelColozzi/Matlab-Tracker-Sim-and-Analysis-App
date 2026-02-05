addpath SignalSources\
addpath Trackers\
addpath Filters\
clear all
threshold = 1;

OnState = NoisedOnState(makedist("Lognormal","mu",0,"sigma",1));
OffState = OffState();

transition = [0.9,0.1;0.01,0.99];
states = {OnState,OffState};

chain = Markov(transition,1,states);

trials = 1000;
data = zeros(trials,1);
currentStates = zeros(trials,1);
for n=1:trials
    chain.Update();
    data(n) = chain.getSignalOutput();
    currentStates(n) = chain.CurrentState;
end

%%
clear all
close all hidden
rng(1234)
N = 100;
trueStates = zeros(4,1,N);
trackStates = zeros(4,N);

trueState = normrnd(0,1,4,1);
dt = 0.1;
A = [[eye(2),dt * eye(2)];[-0.0*eye(2),0.999*eye(2)]];
%Trk = IdealKFTrack(zeros(4,1),eye(4),A,0.01*eye(4),[eye(2),zeros(2)]);
InitialCovPad = zeros(4);
InitialCovPad(3:4,3:4) = 9 * eye(2);
Q = zeros(4);
Q(1:2,1:2) = 0.00001 * eye(2);
Q(3:4,3:4) = 0.01 * eye(2);
Trker = GNN      (0.9,... %PD
                  1e-6,... %Beta %3,... % number of Lags
                  A,Q,[eye(2),zeros(2)],InitialCovPad);

NumTargets = 5;
Targets = cell(1,NumTargets);
Targets = cellfun(@(x)(normrnd(0,1,4,1)),Targets,"UniformOutput",false);
trueStates = horzcat(Targets{:});
for n=1:N
    Measurement = {};
    R = {};
    for targ = 1:NumTargets
    Targets{targ} = A * Targets{targ} + mvnrnd(zeros(4,1),Q,1)';
    Targets{targ}(3:4) = (-1).^(abs(Targets{targ}(1:2))>10).*Targets{targ}(3:4);
    Targets{targ}(1:2) = max(min(Targets{targ}(1:2),10),-10);
    if rand<0.9
        Measurement(end+1) = {Targets{targ}(1:2) + normrnd(0,2,2,1)};
        R(end+1) = {4 * eye(2)};
    end
    end
    %Trk.WeightedUpdate(Measurement,{0.01 * eye(2)},1);
    Trker.UpdateWithMeasurements(Measurement,R);

    trueStates(:,:,n) = horzcat(Targets{:});
end
TruthStates = trueStates(:,:,:);
TrackLog = Trker.LoggingCell;

TrackerApp(TruthStates,TrackLog)

function pd=PD(SNR,T)
    pd = exp(-T./(1+SNR));
end

function pfa = PFA(T)
    pfa = 1-gammainc(T,0);
end