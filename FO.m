%__________________________________________________________________     _%
%  Fertilization optimization algorithm: 
%  New metaheuristic optimization algorithm for  mathematical optimization
%  (FO) source codes version 1.0                                       %
%                                                                        %
%  Developed in MATLAB R2016b                                            %
%                                                                        %
%  Author and programmer: Hazim Nasir Ghafil                             %
%                                                                        %
%         e-Mail: hazimn.bedran@uokufa.edu.iq                            %
%                 hazimbedran@gmail.com                                  %
%                                                                        %
% Homepage: http://staff.uokufa.edu.iq/en/cv.php?hazimn.bedran           %                          
% Main paper:
% Hazim Nasir Ghafil, Shaymaa Alsamia, and Károly Jármai. 
% Fertilization optimization algorithm on CEC2015 and large scale problems
% Pollack Periodica (2021), DOI: https://doi.org/10.1556/606.2021.00343.                     
% https://akjournals.com/view/journals/606/aop/article-10.1556-606.2021.00343/article-10.1556-606.2021.00343.xml%
%___________________________________________________________________     %
% clear;
% clc;
% close all;
%% Problem Definition
% test function
function [sBest,pBest,bestfits] = FO(N,T, lb,ub,dim,fobj)
CostFunction = fobj;  % Objective Function 
nVar = dim;           % Number of Decision Variables
v_size = [1 nVar]; % Decision Variables Matrix Size
% max and min limit
v_min = lb;       % D25
v_max = ub;      % 
%% IWO Parameters
nPop = N;      % Maximum Population Size
MaxIt = T;    % Maximum Number of Iterations

%% Initialization
% Empty Plant Structure
empty_sperm.Position = [];
empty_sperm.velocity = [];
empty_sperm.Cost = [];
sperm = repmat(empty_sperm, nPop, 1);    % Initial Population Array
for i = 1:nPop    
    % Initialize Position
    sperm(i).Position = unifrnd(v_min, v_max, v_size); 
    sperm(i).velocity = sperm(i).Position;
    % Evaluation
    sperm(i).Cost = CostFunction(sperm(i).Position);    
end
% Initialize Best Cost History
bestfits = zeros(MaxIt, 1);
%% IWO Main Loop
x=1;
d=0.95;

% for it = 1:MaxIt  
for it=1:MaxIt
     newsperm = [];    
        for i = 1:nPop
            % Initialize new sperm
           K= (sperm(1).Position+sperm(round(nPop/2)).Position+sperm(end).Position)/3;           
            newsol = empty_sperm;
            L = Levy(nVar);  
            newsol.velocity =  exp(-1/x)*sperm(i).velocity; 
            ds=L.*(sperm(i).Position-sperm(i).velocity);             
            newsol.Position = sperm(i).Position - newsol.velocity+ds-K;           
       
            newsol.Position = max(newsol.Position, v_min);
            newsol.Position = min(newsol.Position, v_max);                
            newsol.Cost = CostFunction(newsol.Position);        
        % Add new sperm to the Population
            newsperm = [newsperm
                      newsol];  %#ok                       
        end
        
    % Merge Populations
    sperm = [sperm
           newsperm ];    
    % Sort Population
    [~, SortOrder]=sort([sperm.Cost]);
    sperm = sperm(SortOrder);

    % Competitive Exclusion (Delete Extra Members)
    if numel(sperm)>nPop
        sperm = sperm(1:nPop);
    end    
    % Store Best Solution Ever Found
    BestSol = sperm(1);    
    % Store Best Cost History
    bestfits(it) = BestSol.Cost;
    x=x*d;
    % Display Iteration Information    

% disp(['It ' num2str(it) ': Best Cost = ' num2str(bestfits(it))]);
end
sBest=BestSol.Cost;
pBest=BestSol.Position;
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end