%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA106
% Project Title: Real-Coded Simulated Annealing in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

% clc;
% clear;
% close all;
function [sBest,pBest,BestCost] = SA(N,T, lb,ub,dim,fobj)
%% Problem Definition

CostFunction=fobj;    % Cost Function

nVar = dim;              % Number of Decision (Unknwon) Variables

VarSize = [1 nVar];     % Decision Variables Matrix Size

VarMin =  lb;           % Lower Bound of Decision Variables
VarMax =  ub;           % Upper Bound of Decision Variables


%% SA Parameters

MaxIt=T;     % Maximum Number of Iterations

MaxSubIt=20;    % Maximum Number of Sub-iterations

T0=0.1;       % Initial Temp.

alpha=0.99;     % Temp. Reduction Rate

nPop=N;        % Population Size

nMove=5;        % Number of Neighbors per Individual

mu = 0.5;       % Mutation Rate

sigma = 0.1*(VarMax-VarMin);    % Mutation Range (Standard Deviation)

%% Initialization

% Create Empty Structure for Individuals
empty_individual.Position=[];
empty_individual.Cost=[];

% Create Population Array
pop=repmat(empty_individual,nPop,1);

% Initialize Best Solution
BestSol.Cost=inf;

% Initialize Population
for i=1:nPop
    
    % Initialize Position
    pop(i).Position=unifrnd(VarMin, VarMax, VarSize);
    
    % Evaluation
    pop(i).Cost=CostFunction(pop(i).Position);
    
    % Update Best Solution
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
    
end

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Intialize Temp.
T=T0;

%% SA Main Loop

for it=1:MaxIt
    
    for subit=1:MaxSubIt
        
        % Create and Evaluate New Solutions
        newpop=repmat(empty_individual,nPop,nMove);
        for i=1:nPop
            for j=1:nMove
                
                % Create Neighbor
                newpop(i,j).Position=Mutate(pop(i).Position,mu,sigma,VarMin,VarMax);
                
                % Evaluation
                newpop(i,j).Cost=CostFunction(newpop(i,j).Position);
                
            end
        end
        newpop=newpop(:);
        
        % Sort Neighbors
        [~, SortOrder]=sort([newpop.Cost]);
        newpop=newpop(SortOrder);
        
        for i=1:nPop
            
            if newpop(i).Cost<=pop(i).Cost
                pop(i)=newpop(i);
                
            else
                DELTA=(newpop(i).Cost-pop(i).Cost)/pop(i).Cost;
                P=exp(-DELTA/T);
                if rand<=P
                    pop(i)=newpop(i);
                end
            end
            
            % Update Best Solution Ever Found
            if pop(i).Cost<=BestSol.Cost
                BestSol=pop(i);
            end
        
        end

    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
%     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    % Update Temp.
    T=alpha*T;
    
    sigma = 0.98*sigma;
    
end

sBest=BestSol.Cost;
pBest=BestSol.Position;

end

function y=Mutate(x,mu,sigma,VarMin,VarMax)

    A=(rand(size(x))<=mu);
    J=find(A==1);

    y=x;
    y(J)=x(J)+sigma*randn(size(J));

    % Clipping
    y=max(y,VarMin);
    y=min(y,VarMax);
    
end