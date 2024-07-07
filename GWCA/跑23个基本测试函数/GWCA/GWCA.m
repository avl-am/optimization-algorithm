%%
function [Worker1_fit,Worker1,Convergence_curve]=GWCA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
    neval = 0;
    Worker1=zeros(1,dim);   
    Worker1_fit=inf;
    Worker2=zeros(1,dim);   
    Worker2_fit=inf;
    Worker3=zeros(1,dim);   
    Worker3_fit=inf;
    %Initializing
    SL=1;
    T=8.3;
    g=9.8;
    m=3;
    e=0.1;
    P=9;
    Q=6;
    Cmax=exp(3);
    Cmin=exp(2);
    N=30;
    %Initialize the positions of search agents
    Positions=initial(N,dim,ub,lb);
    LNP=ceil(N*e);
    Convergence_curve=zeros(1,Max_iter);
    
    l=0;% Loop counter
    % Main loop
    % Max_iter = Max_iter
    while l < Max_iter
        %% Initialize the population position and the best fitness value
        Positions=initial(SearchAgents_no,dim,ub,lb);                           %Eq.(1)
        for i=1:N
            Fitness(i,:)=fobj(Positions(i,:));
            if Fitness(i)<Worker1_fit
                Worker1_fit=Fitness(i);  Worker1=Positions(i,:);
            elseif Fitness(i)>Worker1_fit && Fitness(i)<Worker2_fit
                Worker2_fit=Fitness(i);  Worker2=Positions(i,:);
            elseif Fitness(i)>Worker1_fit && Fitness(i)>Worker2_fit && Fitness(i)<Worker3_fit
                Worker3_fit=Fitness(i);  Worker3=Positions(i,:);
            end
        end
        Pbest=Fitness;
        P_Pobest=Positions;
        
        %% Start
        for it=1:size(Positions,1)  
            for i=1:size(Positions,2)  
                Index=randi([1,3]);
                if Index==1
                    %% engineer movement
                    C=log((Cmax-Cmin)*(Max_iter-it)/Max_iter+Cmin);              %Eq.(7)
                    H=(1-it/Max_iter);                                          %Eq.(6)
                    sitar=80*rand(1,dim);
                    TL=SL-it/Max_iter+eps;
                    a=(T.*TL)./m-g.*(H./sin(sitar));
                    v=a.*C.*gampdf(it,P,Q);                                    %Eq.(5)
                    Study=(-1)^randi([0,1])*(Worker1-Positions(i,:)).*rand(1,dim);%Eq.(4)
                    Positions(i,:)=Worker1+Study+Positions(i,:).*v.*rand(1,dim);                            %Eq.(3)
                elseif Index==2
                    %% soldier movement
                    C=log((Cmax-Cmin)*(Max_iter-it)/Max_iter+Cmin);
                    H=(1-it/Max_iter)+eps;
                    sitar=80*rand(1,dim);
                    v=m.*g.*(H./sin(sitar)).*C.*gampdf(it/Max_iter,P,Q);        %Eq.(9)
                    Index=1:N;
                    Index(Index==i)=[];
                    a=Fitness(Index)-Fitness(i);
                    [~,ide]=min(abs(a));
                    improve=sign(Fitness(ide)-Fitness(i))*(Positions(ide,:)-Positions(i,:)).*v.*rand(1,dim);
                    study=(Worker2-Positions(i,:)).*rand(1,dim);
                    Positions(i,:)=Positions(i,:)+improve+study;                             %Eq.(10)
                elseif Index==3
                    %% labor movement
                    Positions(i,:)=Positions(i,:)+2*(Worker3-Positions(i,:)).*rand(1,dim)+(P_Pobest(i,:)-Positions(i,:)).*gampdf(it/Max_iter,P,Q);                   %Eq.(11)
                end
                %% constraint bound
                Flag4ub=Positions(i,:)>ub;
                Flag4lb=Positions(i,:)<lb;
                Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
                Fitness(i,:)=fobj(Positions(i,:));
                neval = neval + 1;
                %update local optimum
                if Fitness(i,:)<Pbest(i,:)
                    Pbest(i,:)=Fitness(i,:);
                    P_Pobest(i,:)=Positions(i,:);
                end
                %update leader
                if Fitness(i)<Worker1_fit
                    Worker1_fit=Fitness(i);  
                    Worker1=Positions(i,:);

                elseif Fitness(i)>Worker1_fit && Fitness(i)<Worker2_fit
                    Worker2_fit=Fitness(i);  
                    Worker2=Positions(i,:);

                elseif Fitness(i)>Worker1_fit && Fitness(i)>Worker2_fit && Fitness(i)<Worker3_fit
                    Worker3_fit=Fitness(i);  Worker3=Positions(i,:);
                end
            end
            %% Personnel elimination mechanism
            [~,Index]=sort(Fitness,'descend');
            Positions(Index(1:LNP),:)=(ub-lb).*rand(LNP,dim)+lb;
            %% sort fitness
            record_fitness(it)=Worker1_fit;
            disp(['Iter',num2str(it),'   Best Fitness',num2str(Worker1_fit)])
        end
        %% Record the fitness value curve and optimal parameters obtained each time, as well as the best fitness value
        % Nfe(irun,:) = neval;
        % Conver(irun,:)=record_fitness;
        % DestinationFitness(irun,:)=Worker1_fit;
        % BestPosition(irun,:)=Worker1;
        l=l+1;    
        Convergence_curve(l)=Worker1_fit;
    end
    % BestData.Congervence = mean(record_fitness,1);
    % BestData.MeanNfe = mean(Nfe);
    % BestData.Conver=mean(Conver,1);
    % BestData.ALLFitness=DestinationFitness;
    % BestData.StdScore=std(DestinationFitness);
    % BestData.MeanScore=mean(DestinationFitness);
    % [BestData.Fitness,MINindex]=min(DestinationFitness);
    % BestData.BestPosition=BestPosition(MINindex,:);
end

function Y=gampdf(x,a,b)
    Y = 1/(gamma(x)*b^a)*x^(a-1)*exp(-x/b);
end


%% Level 1: Initializing
function Positions=initial(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries
cxl=rand(SearchAgents_no,dim);
for j=1:dim
    if cxl(j)==0
        cxl(j)=0.1;
    end
    if cxl(j)==0.25
        cxl(j)=0.26;
    end
    if cxl(j)==0.5
        cxl(j)=0.51;
    end
    if cxl(j)==0.75
        cxl(j)=0.76;
    end
    if cxl(j)==1
        cxl(j)=0.9;
    end
end
for j=1:dim
    cxl(j)=4*cxl(j)*(1-cxl(j));        %logic混沌方程
end
% if Boundary_no==1
    Positions=cxl.*(ub-lb)+lb;
% end
% If each variable has a different lb and ub
% if Boundary_no>1
%     for i=1:dim
%         ub_i=ub(i);
%         lb_i=lb(i);
%         Positions(i,:)=cxl.*(ub_i-lb_i)+lb_i;
%     end
end


