classdef DefaultObserver_CableLengthNoise < handle
    %DEFAULTOBSERVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Xhat_p
        Xhat_pDOT
        Xhat_RL
        Xhat_L
    end
    
    methods
        function obj = DefaultObserver_CableLengthNoise(~,~,~,~,~,~)
            %do nothing
        end
        
        function [Xhat,Z,Pm,Pp] = estimateState(obj,Y,~,~,~)
            %ESTIMATESTATE This function takes in the full-state
            %observations from 'defaultPlant' and parses the struct for the
            %necessary Xhat,Uhat output for the controller
            
            %Example of recording down current observations as a property,
            %maybe useful for observers that require updating based on past
            %estimations (e.g., Kalman Filter)
            obj.Xhat_p = Y.p;
            obj.Xhat_pDOT = Y.pDOT;
            obj.Xhat_RL = Y.RL;
            obj.Xhat_L = Y.L;
            
            Xhat.p = obj.Xhat_p;
            Xhat.pDOT = obj.Xhat_pDOT;
            Xhat.RL = obj.Xhat_RL;
            Xhat.L = obj.Xhat_L;
            Z = zeros(84,1);
            Pm = zeros(numel(Y.p)+numel(Y.pDOT)+numel(Y.RL)+numel(Y.L),1);
            Pp = zeros(numel(Y.p)+numel(Y.pDOT)+numel(Y.RL)+numel(Y.L),1);
            
            %inject noise into cable lengths
            frac = 0.2; %range from min to max percent deviation
            Xhat.RL = Xhat.RL.*(1+frac*(rand(size(Xhat.RL))-0.5));
        end
    end
end

