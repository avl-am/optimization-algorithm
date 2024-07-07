function [fbest, xbest, Convergence_curve] = LPO(N,maxEvals, lb, ub, dim,fobj)

    nVar=dim;                     % Number of Decision Variables
    VarSize = [1 nVar];          % Decision Variables Matrix Size
    VarMin = lb;               % Lower Bound of Decision Variables
    VarMax = ub;            % Upper Bound of Decision Variables

    %% LPO Parameters
    MaxIt = maxEvals;              % Maximum Number of Iterations
    nPop =N;                    % Maximum Population Size

    %% Initialization
    % Empty Plant Structure
    empty_plant.Position = [];
    empty_plant.Cost = [];
    pop = repmat(empty_plant, nPop, 1);    % Initial Population Array
    BestSol.Cost=inf;

    for i = 1:numel(pop)
        % Initialize Position
        pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
        Delta(i)=rand*2*pi;
        pop(i).Cost=(fobj(pop(i).Position));%if you want to use form  Cost Functions CEC
        sigma1(i,:) =rand(VarSize);
        if pop(i).Cost<=BestSol.Cost
            BestSol.Position=pop(i).Position;
            BestSol.Cost=pop(i).Cost;
        end
        BestCost1(i)=BestSol.Cost;
    end

    it=nPop;
    % % for it=1:MaxIt
    for it = 1:MaxIt

        % Get Best and Worst Cost Values
        Costs = [pop.Cost];
        BestCost = min(Costs);
        WorstCost = max(Costs);

        newpop = [];

        for i = 1:numel(pop)

            R=[];
            C=[];
            for jj = 1:5

                newsol = empty_plant;
                newsol2 = empty_plant;
                R= pop(i).Cost;               %% Eq(9)
                C=(R/2)*sin(Delta(i));        %% Eq(10)

                
                if jj==1 %% Eq(7)
                    newsol.Position =pop(i).Position+ (((R^2+(1/(2*pi*nVar*R*C)^2))^-0.5)*sin(2*pi*nVar*it)*sin((2*pi*nVar*it)+Delta(i))).*pop(i).Position;
                else     %% Eq(16)
                    newsol.Position =pop(i).Position+ (((R^2+(1/(2*pi*nVar*R*C)^2))^-0.5)*sin(2*pi*nVar*it)*sin((2*pi*nVar*it)+Delta(i))).*(Pos1);
                end

                %%Eqs(12)-(13)
                A1=randperm(nPop);
                A1(A1==i)=[];
                a1=A1(1);
                a2=A1(2);
                a3=A1(3);

                Costs = [pop.Cost];
                BestCost = min(Costs);
                WorstCost = max(Costs);

                %%Eq(13)
                aa1=((pop(a2).Cost-pop(a3).Cost)/abs(pop(a3).Cost-pop(a2).Cost));
                if (pop(a2).Cost-pop(a3).Cost)==0
                    aa1=1;
                end
                %%Eq(13)
                aa2=((pop(a1).Cost-pop(i).Cost)/abs(pop(a1).Cost-pop(i).Cost));
                if (pop(a1).Cost-pop(i).Cost)==0
                    aa2=1;
                end
                
                %%Mnew,2 in Eq(12)
                newsol.Position= newsol.Position+(aa2*(sigma1(i))).*(newsol.Position-pop(a1).Position)+(aa1*sigma1(i)).*(pop(a3).Position-pop(a2).Position);
                
                %%Mnew,3 in Eq(14)
                newsol2.Position =pop(a1).Position+(sigma1(i)).*(pop(a3).Position-pop(a2).Position);
                see=pop(i).Position;

                %%Eq(14)
                for j =1:nVar
                    if rand/jj>rand  %% Eq(15)
                        Pos1(1,j) =newsol2.Position(1,j);
                    else
                        Pos1(1,j)=newsol.Position(1,j);
                    end
                end


                % Apply Lower/Upper Bounds
                Pos1 = max(Pos1, VarMin);
                Pos1 = min(Pos1, VarMax);
                newsol.Position=Pos1;
                Delta(i)=atan((1/(2*pi*nVar*R*C)));  % Eq(8)

                newsol.Cost=(fobj(newsol.Position));%if you want use form  Cost Functions CEC 2014
                if newsol.Cost<pop(i).Cost
                    pop(i)=newsol;
                    if pop(i).Cost<=BestSol.Cost
                        BestSol=pop(i);
                    end
                end


                BestCost1(it)=BestSol.Cost;%if you want use form  Cost Functions CEC 2014
                sigma1(i,:) =rand(VarSize);
            end
        end

        %Store Best Cost History
        Convergence_curve(it) = BestSol.Cost;
        fbest = BestSol.Cost;
        xbest = BestSol.Position;

        % %Display Iteration Information
        % if mod(it,180)==0
        %     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost)]);
        % end
    end


    % %% Results
    % Cost_Rsult(1, jin)=BestSol.Cost;%if you want use form  Cost Functions CEC 2014 or CEC2017
    % Rsult(jin,:)=BestCost1;
end






