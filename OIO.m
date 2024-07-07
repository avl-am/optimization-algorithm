% % % % This codes is related to the Optics Inspired Optiization and its variant for the following problem
% % % % min f(X)
% % % % Xmin<=X<=Xmax
% % % % You just need to introduce your function to be minimized in objfn.m. Then enter the problem parameters
% % % % (lines 20 to 23) and OIO parameters(lines 26 to 40), accordingly. OIO parameters need to be tuned based
% % % % on the problem on hand. This code is run under atmost two subpopulations. You can choose to have one population
% % % % or two subpopulations. The number of changes in the new image solution is controlled with "NumberofChanges" for
% % % % the one population case and is controlled with "NumberofChangesforPop1" and "NumberofChangesforPop2" for the
% % % % two subpopulation case. When NumberofChanges=N or NumberofChangesforPop2=N, it is recomended to use a
% % % % perturbation term in the OIO equation which is activated under is_random_perturbation_required=1
function [ BestCost,BestSolution,best_objective] = OIO(npop,maxit,Xmin,Xmax,N,objfun)

% global CostFunction VarSize VarMin VarMax MaxIt nPop NFE nvar AlgorithmTermination

% format long
% clear all
best_objective=[];
num_of_function_evaluations=[];
execution_time=[];
%%%%%%%%%%%%%%%%%%%problem parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target_value_to_reach=inf;
if length(Xmin)==1
    Xmin=Xmin*ones(1,N);  %%a lower bound value on the optimization parameters
    Xmax=Xmax*ones(1,N);  %%an upper bound value on the optimization parameters
end
%%%%%%%%%%%%%%%%%%%%inputs parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumofLightPoints=npop;%should be devided into 3 and greater than or equal to 9
type_of_algorithm=1;      % 1 for OIO         2 for ROIO        3 for COIO
NumberofPopulations=2;    % 1  or 2
if NumberofPopulations==1
    NumberofChanges=N;   %% an integer between 1 and N
    SizeofPopulation1=NumofLightPoints;
elseif NumberofPopulations==2
    SizeofPopulation1=floor(NumofLightPoints/3);
    NumberofChangesforPop1=1; %% an integer between 1 and N
    NumberofChangesforPop2=N; %% an integer between 1 and N
end
is_random_perturbation_required=1;   %% 0 is no 1 is yes. 1 is recomended when NumberofChanges=N or  NumberofChangesforPop2=N
NumberOfMirrors_min=1;
NumberOfMirrors_max=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=cputime;
stopflag=0;
continueflag=0;
eval=0;
globalbest=realmax*ones(1,N+1);
%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%
initialpop=[];
for pop=1:N
    initialpop(:,pop)=unifrnd(Xmin(pop),Xmax(pop),NumofLightPoints,1);
end
cost=[];
for i=1:NumofLightPoints
    LightPoint=initialpop(i,:);
    objectivefunction=objfun(LightPoint);
    
    cost=[cost
        objectivefunction];
    if objectivefunction<globalbest(end)
        globalbest=[LightPoint objectivefunction];
    end
end
initialpopulation=[initialpop cost];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binahayat=abs(max(initialpopulation(:,end)));
BestLightPointPosition=initialpopulation;
iteration=0;
while iteration<maxit
    for counter=1:NumofLightPoints
        currentlightpoint=BestLightPointPosition(counter,:);
        componentMat=[];
        NumberOfMirrors=randi([NumberOfMirrors_min,NumberOfMirrors_max],1,1);
        for Mirror=1:NumberOfMirrors
            the_mirror_is_concave=0;
            the_mirror_is_convex=0;
            %%%%%%%%%%%%%  roulete wheel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            difendmat=[[1:NumofLightPoints]' repmat(currentlightpoint(end),NumofLightPoints,1)-BestLightPointPosition(:,end)];
            candidsol=find(difendmat(:,end)~=0);
            gbid=min(find(BestLightPointPosition(:,end)==min(BestLightPointPosition(:,end))));
            if counter<=SizeofPopulation1
                candidsol=candidsol(find(candidsol<=SizeofPopulation1));
            else
                candidsol=candidsol(find(candidsol>SizeofPopulation1));
            end
            if size(find(candidsol==gbid),1)==0
                if currentlightpoint(end)==min(BestLightPointPosition(:,end))
                else
                    candidsol=[candidsol
                        gbid];
                end
            end
            if size(candidsol,1)==0
                if NumberofPopulations==1 || sum(repmat(max(BestLightPointPosition(:,end)),NumofLightPoints,1)-BestLightPointPosition(:,end))==0
                    stopflag=1;
                    break
                else
                    continueflag=1;
                    continue
                end
            end
            roulettecandid=BestLightPointPosition(candidsol,end);
            pos=find(roulettecandid>=0);
            neg=find(roulettecandid<0);
            roulettecandid(pos)=1./(1+roulettecandid(pos));
            roulettecandid(neg)=1+abs(roulettecandid(neg));
            refinedvalues=roulettecandid;
            probabilities=(refinedvalues./sum(refinedvalues));
            select=cumsum(probabilities)-repmat(rand,size(roulettecandid,1),1)>0;
            focus=candidsol(min(find(select==1)),1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XO=currentlightpoint;
            FO=BestLightPointPosition(focus,:);
            if XO(end)>FO(end)
                the_mirror_is_concave=1;
                F=[BestLightPointPosition(focus,1:N) unifrnd(XO(end),XO(end)+binahayat)];
                X=[currentlightpoint(1:N) unifrnd(XO(end),XO(end)+binahayat)];
            else
                the_mirror_is_convex=1;
                F=[BestLightPointPosition(focus,1:N) unifrnd(XO(end)-binahayat,XO(end))];
                X=[currentlightpoint(1:N) unifrnd(FO(end),FO(end)+binahayat)];
            end
            %%spherical abberation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
            R=(F(end)-FO(end));
            HO=norm(X(1:end-1)-F(1:end-1));
            while (HO>abs(R)) || (((abs(R)^2)/(2*sqrt((abs(R)^2)-(HO^2))))-(abs(R)/2)>0.01)
                binahayat=2*binahayat;
                if the_mirror_is_concave==1
                    F(end)=F(end)+binahayat;
                elseif the_mirror_is_convex==1
                    F(end)=F(end)-binahayat;
                end
                R=(F(end)-FO(end));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of sperical aberration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p=(X(end)-FO(end));
            q=(p*R)/(2*p-R);
            Hi=-q*HO/p;
            if type_of_algorithm==1
                componentMat=[componentMat
                    F(1:end-1)-(q/p)*(X(1:end-1)-F(1:end-1))];
            elseif type_of_algorithm==2
                ref=(2*F(1:end-1)*X(1:end-1)'/(norm(X(1:end-1))^2))*X(1:end-1)-F(1:end-1);
                componentMat=[componentMat
                    F(1:end-1)-HO*(q/p)*((F(1:end-1)-ref)/norm((F(1:end-1)-ref)))];
            elseif type_of_algorithm==3
                ref=(2*F(1:end-1)*X(1:end-1)'/(norm(X(1:end-1))^2))*X(1:end-1)-F(1:end-1);
                convex=rand;
                componentMat=[componentMat
                    F(1:end-1)-HO*(q/p)*(convex*((X(1:end-1)-F(1:end-1))/HO)+(1-convex)*((F(1:end-1)-ref)/norm((F(1:end-1)-ref))))];
            end
            if is_random_perturbation_required==1
                if NumberofPopulations==1
                    componentMat(end,:)=componentMat(end,:)+HO*(rand-1/2)*binornd(1,1-(iteration/maxit));
                elseif NumberofPopulations==2  &&  counter>SizeofPopulation1
                    componentMat(end,:)=componentMat(end,:)+HO*(rand-1/2)*binornd(1,1-(iteration/maxit));
                end
            end
        end
        
        if continueflag==1;
            continueflag=0;
            continue
        end
        if NumberOfMirrors==1
            newprojection=componentMat;
        elseif NumberOfMirrors==2
            Randomcoefficients=rand;
            newprojection=Randomcoefficients*componentMat(1,:)+(1-Randomcoefficients)*componentMat(2,:);
        else
            Randomcoefficients=rand;
            for ran=1:NumberOfMirrors-2
                Randomcoefficients=[Randomcoefficients unifrnd(0,1-sum(Randomcoefficients))];
            end
            Randomcoefficients=[Randomcoefficients 1-sum(Randomcoefficients)]';
            newprojection=sum(repmat(Randomcoefficients,1,N).*componentMat);
        end
        outofrangemin=find(newprojection(1:N)<Xmin);
        outofrangemax=find(newprojection(1:N)>Xmax);
        newprojection(1,outofrangemin)=Xmin(1,outofrangemin)+rand(1,size(outofrangemin,2)).*(Xmax(1,outofrangemin)-Xmin(1,outofrangemin));
        newprojection(1,outofrangemax)=Xmin(1,outofrangemax)+rand(1,size(outofrangemax,2)).*(Xmax(1,outofrangemax)-Xmin(1,outofrangemax));
        SELECT=randperm(N);
        if NumberofPopulations==2
            if counter<=SizeofPopulation1
                NumberofChanges=NumberofChangesforPop1;
            else
                NumberofChanges=NumberofChangesforPop2;
            end
        end
        CD=SELECT(1,1:NumberofChanges);
        changed_projection=currentlightpoint(1,1:N);
        changed_projection(1,CD)=newprojection(1,CD);
        LightPoint=changed_projection;
        objectivefunction=objfun(LightPoint);
        eval=eval+1;
        if objectivefunction<globalbest(end)
            globalbest=[changed_projection  objectivefunction];
            %                  fprintf('BestObj=%d eval=%g\n',objectivefunction,eval);
        end
        
        if objectivefunction<BestLightPointPosition(counter,end)
            BestLightPointPosition(counter,:)= [changed_projection  objectivefunction];
        end
    end
    iteration=iteration+1;
    best_objective (iteration) = globalbest(end);
   
    
end
BestCost=globalbest(end);
BestSolution=globalbest(1:end-1);

end
