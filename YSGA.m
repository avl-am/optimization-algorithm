
%% Yellow Saddle Goatfish Algorithm (YSGA)                                    
%       Source Paper:
%           Daniel Zald¨ªvar, Bernardo Morales, Alma Rodr¨ªguez, Arturo Valdivia-G, Erik Cuevas and Marco P¨¦rez-Cisneros
%           "A novel bio-inspired optimization model based on Yellow Saddle Goatfish behavior"
%           https://doi.org/10.1016/j.biosystems.2018.09.007
%
%       Coded by: Bernardo Morales
%       Questions about the code: jb.moralescastaneda@gmail.com 
%%  FUNCTION INSTRUCTIONS:
%   This function allows to implement the YSGA algorithm for the
%   optimization of an user defined function. All other required
%   sub-functions are included in this script. To implement YSGA,
%   use the following command line as a template:
%
%   [bestFit,bestFitHist,bestSol] = YSGA(func,population,dims,itern, lb, ub)
%   
%   INPUTS:
%       population   -> Population Size
%       itern        -> Maximum number of iterations 
%       func         -> Objective Function. Invoke Function Handle for this (i.e.@Ackley)
%       dims         -> Search Space total dimensions
%       lb           -> Search Space Lower Bound 
%       ub           -> Search Space Upper Bound      
%   OUTPUTS:
%       bestSol       -> Best Found Solution
%       bestFit       -> Fitness (quality) of the best found solution
%       bestFitHist   -> Best Solutions Fitness Evolution by iterations
%
function [bestFit,bestSol,bestFitHist] = YSGA(population,itern, lb, ub,dims,func)
% Initialization
n = population;  
maxGeneration = itern;
dim = dims;
groupN = 4;
range = [lb, ub];
f = func;
globalBest.Sol = zeros(1,dim);
globalBest.Fitness = Inf;
globalBest.Group = 0;
fish.Pos = cell(1,groupN);
fish.Fitness = cell(1,groupN);
fish.groupBest = zeros(1,groupN);
FitHist = zeros(1,maxGeneration);
flags = zeros(1,groupN);
[fish] = ysga_init(f,fish,n,range,groupN,dim);
% Main iterative process
for i = 1:maxGeneration
    % For each group of fish
    for j=1:groupN
       zn = cell2mat(fish.Fitness(1,j));
       alpha = -1+i*((-1)/maxGeneration);
       % Chaser fish movement
       [fish] = ysga_hunt(fish,j,range,f,dim,i,maxGeneration,globalBest);
       % Blocker fish movement
       [fish] = ysga_block(fish,j,range,alpha,dim);  
       xn2 = cell2mat(fish.Pos(1,j));
       z2 = zeros(size(xn2(:,1)));
       ni = size(z2);
       for k = 1:ni(1) 
           z2(k) = f(xn2(k,:)); 
       end
       fish.Fitness{1,j} = z2;
       [zn2,Ind2] = sort(cell2mat(fish.Fitness(1,j)), 'ascend');
       for k=1:dim
            xn2(:,k) = xn2(Ind2,k);
       end
       fish.Pos{1,j} = xn2;
       fish.Fitness{1,j} = zn2;
       %Update best found global solution
       if zn2(1) < globalBest.Fitness
           globalBest.Fitness = zn2(1);
           globalBest.Sol = xn2(1,:);
           globalBest.Group = j;
       %Update flag for change of zone    
       elseif zn2(1) > globalBest.Fitness
           flags(j) = flags(j) + 1; 
       end
       %Update best found solution of the current group
       if zn2(1) < zn(1)
           fish.groupBest = zn2(1);
       elseif zn2(1)== zn(1)
           flags(j) = flags(j) + 1; 
       end
    end
    
    for j=1:groupN
        if flags(j) >= 10 && j ~= globalBest.Group
            %Start change zone behaviour
            ysga_changezone(fish,j,range,dim,globalBest);
            flags(j) = 0;
        end
    end
   FitHist(i) = globalBest.Fitness;
end  
bestFit = globalBest.Fitness;
bestSol = globalBest.Sol;
bestFitHist = FitHist;
function [fish] = ysga_init(f,fish,n,range,groupN,dim)
xo = zeros(dim,n);
for i=1:dim
    xo(i,:) = rand(1,n)*(range(2)-range(1))+range(1);
end
[idx,~] = kmeans(transpose(xo),groupN,'Distance','sqeuclidean','Replicates',5);
for j=1:groupN
    fish.Pos{1,j} = transpose(xo(:,idx==j));
    xo2 = cell2mat(fish.Pos(1,j));
    zo2 = zeros(size(xo2(:,1)));
    ni = size(zo2);
    for k = 1:ni(1) 
        zo2(k) = f(xo2(k,:)); 
    end
    fish.Fitness{1,j} = zo2;
    [zn,Ind]=sort(cell2mat(fish.Fitness(1,j)), 'ascend');
    xn = zeros(size(xo2));
    for k=1:dim
         xn(:,k) = xo2(Ind,k);
    end
    fish.Pos{1,j} = xn;
    fish.Fitness{1,j} = zn;
end
function [fish] = ysga_hunt(fish,j,range,f,dim,iter,maxGeneration,globalBest)
xo = transpose(cell2mat(fish.Pos(1,j)));
zo = transpose(cell2mat(fish.Fitness(1,j)));
xn = zeros(dim,10);
zn = zeros(1,10);
beta = 1.99 + (.001*iter/(maxGeneration/10));
if beta > 2
    beta = 2;
end
levySteps = levyStep(10,dim,beta);
 if (transpose(xo(:,1))-globalBest.Sol) ~= 0
    for i = 1:10
         levySteps(i,:) = 1*levySteps(i,:).*(transpose(xo(:,1))-globalBest.Sol);
         xn(:,i) = xo(:,1) + transpose(levySteps(i,:));  
    end
 else
    for i = 1:dim
        xn(i,:) = xo(i,1) + transpose(levySteps(:,i));  
    end
  end
         
for i = 1:10
   zn(i) = f(transpose(xn(:,i)));
end
[Fitness,Ind]=sort(zn, 'ascend');
xn = xn(:,Ind);
zn = Fitness;
if zn(1) < zo(1)
    for i = 1:dim
        xo(i,1) = xn(i,1);
    end
    zo(1) = zn(1);
end
[xo] = findrange(xo,range,dim);
fish.Pos{1,j} = transpose(xo);
fish.Fitness{1,j} = zo;
function [fish] = ysga_block(fish,j,range,alpha,dim)
xo = transpose(cell2mat(fish.Pos(1,j)));
xn = xo;
ni = size(xn(1,:),2);
b = 1;
for i = 2:ni
    t = (alpha-1)*rand+1;
    rx = zeros(1,dim);
    for k = 1:dim
       r = (rand*2-1);
       rx(k) = abs(xn(k,1) * r - xn(k,i)); 
       xn(k,i) = rx(k) * exp(b.*t) .* cos(t.*2*pi) + xn(k,1);
    end
end
[xn] = findrange(xn,range,dim);
fish.Pos{1,j} = transpose(xn);
function [fish] = ysga_changezone(fish,j,range,dim,globalBest)
xo = transpose(cell2mat(fish.Pos(1,j)));
for i=1:dim
     xo(i,:) =  (globalBest.Sol(i) + xo(i,:))/2;
end
[xo] = findrange(xo,range,dim);
fish.Pos{1,j} = transpose(xo);
function [xo] = findrange(xo,range,dim)
for i=1:dim
    for j=1:length(xo(1,:))
        if xo(i,j)<=range(1)
            xo(i,j)=range(1);
        end
        if xo(i,j)>=range(2)
            xo(i,j)=range(2);
        end
    end
end
function [z] = levyStep(n,m,beta)
% This function implements Levy's flight. 
% For more information see 
%'Multiobjective cuckoo search for design optimization Xin-She Yang, Suash Deb'. 
% Coded by Hemanth Manjunatha on Nov 13 2015.
% Input parameters
% n     -> Number of steps 
% m     -> Number of Dimensions 
% beta  -> Power law index  % Note: 1 < beta < 2
% Output 
% z     -> 'n' levy steps in 'm' dimension
    num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator 
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator
    sigma_u = (num/den)^(1/beta);% Standard deviation
    u = random('Normal',0,sigma_u^2,n,m); 
    
    v = random('Normal',0,1,n,m);
    z = u./(abs(v).^(1/beta));
    
     
% Copyright (c) 2016, Hemanth
% All rights reserved.


