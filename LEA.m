
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Love Evolution Algorithm (LEA) source codes version 1.0
%  
%  Developed in:	MATLAB 9.13 (R2022b)
%  
%  Programmer:		Yuansheng Gao (e-mail: gaoyuansheng2021@163.com)
%  
%  Original paper:	Yuansheng Gao, Jiahui Zhang, Yulin Wang
%					Jinpeng Wang, Lang Qin.
%					Love Evolution Algorithm: a stimulus-value-role theory 
%                   inspired evolutionary algorithm for global optimization, 
%                   The Journal of Supercomputing.
%             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H_G,G, CV] = LEA (N,MaxFEs,lb,ub,dim,fun)
%% Initialization
CV  = zeros (1, MaxFEs);                 % Convergence curves
[X,lb,ub] = Initialization(N,dim,ub,lb); % Initialize the population
FE = 0;                                  % Number of function evaluations
H =  zeros (N, 1);                       % The happiness degrees
H_G = inf;                               % The best value

for i = 1:N
    H(i) = fun (X(i,:));
    FE = FE+1;
    if H_G>H(i)
        H_G = H(i,:);
        G = X(i,:);                      % The best solution
    end
    CV(FE) = H_G;
end
% Parameter settings
h_max = 0.7;
h_min = 0;
lambda_c = 0.5;
lambda_p = 0.5;
%% Main loop
while FE<MaxFEs
    h = (1-FE/MaxFEs)*(h_max-h_min)+h_min; %Eq. (17)
    %% Encounter
    % Eqs. (4) and (5)
    r = randperm(N);
    A = X(r(1:N/2),:);
    B = X(r(N/2+1:end),:);
    H_A = H(r(1:N/2));
    H_B = H(r(N/2+1:end));
    %% Stimulus phase
    c = Gap_P(H_A,H_B); % Eq. (6)
    mu = sum(sqrt(sum((X-G).^2)/N))/dim + eps; % Eq. (10)
    for i = 1:N/2
        if c(i)<lambda_c
            %% Value phase
            for j = 1:dim
                % Eq. (11)/(12)
                phi1 = G(:,j)*A(i,j);
                phi2 = G(:,j)^2 + A(i,j)*B(i,j);
                phi3 = G(:,j)*B(i,j);
                % Eq. (13)
                rho_A = sqrt((phi2-phi1)^2);
                rho_B = sqrt((phi2-phi3)^2);
                % Eq. (14)
                A(i,j) = rand* A(i,j) + randn*rho_A;
                B(i,j) = rand* B(i,j) + randn*rho_B;
            end
            % Eq. (19)
            FE = FE+1;if FE>MaxFEs;break;end
            [A(i,:),H_A(i),G,H_G,CV]=Update_A_mod(A(i,:),CV,FE,G,H_G,ub,lb,fun);
            FE = FE+1;if FE>MaxFEs;break;end
            [B(i,:),H_B(i),G,H_G,CV]=Update_B_mod(B(i,:),CV,FE,G,H_G,ub,lb,fun);

            p(i) = (rand+0.5)*c(i)*sum(sqrt((A(i,:)-B(i,:)).^2))/(dim*mu); % Eq. (15)
            if p(i)<lambda_p
                %% Role phase
                % Eq. (16)
                xi = A(i,:).*B(i,:);
                xi = (xi-min(xi))/(max(xi)-min(xi) + eps) + h;
                for j = 1:dim
                    % Eq. (18)
                    A(i,j) = G(:,j) + randn*mu*xi(j);
                    B(i,j) = G(:,j) + randn*mu*xi(j);
                end
            else
                %% Reflection operation
                for j = 1:dim
                    % Eq. (7)
                    sA = (3*rand-1.5)*(A(i,j)/(B(i,j) + eps));
                    sB = (3*rand-1.5)*(B(i,j)/(A(i,j) + eps));
                    % Eq. (8)
                    z = randi(dim);
                    k = randi(dim);
                    delta = 0.5*(A(i,z)/(ub(z)-lb(z)) + B(i,k)/(ub(k)-lb(k)));
                    % Eq. (9)
                    A(i,j) = G(:,j) + sA*mu*delta;
                    B(i,j) = G(:,j) + sB*mu*delta;
                end
            end
        else
            %% Reflection operation
            for j = 1:dim
                % Eq. (7)
                sA = (3*rand-1.5)*(A(i,j)/(B(i,j) + eps));
                sB = (3*rand-1.5)*(B(i,j)/(A(i,j) + eps));
                % Eq. (8)
                z = randi(dim);
                k = randi(dim);
                % Eq. (9)
                delta = 0.5*(A(i,z)/(ub(z)-lb(z)) + B(i,k)/(ub(k)-lb(k)));
                A(i,j) = G(:,j) + sA*mu*delta;
                B(i,j) = G(:,j) + sB*mu*delta;
            end
        end
        % Eq. (20)
        FE = FE+1;if FE>MaxFEs;break;end
        [A(i,:),H_A(i),G,H_G,CV]=Update_A_ordinary(A(i,:),CV,FE,G,H_G,ub,lb,fun);
        FE = FE+1;if FE>MaxFEs;break;end
        [B(i,:),H_B(i),G,H_G,CV]=Update_B_ordinary(B(i,:),CV,FE,G,H_G,ub,lb,fun);
    end
    X = [A;B];
    H = [H_A;H_B];
end
end

% Eq. (6)
function p = Gap_P(f1,f2)
p = (0.5+rand(length(f1),1)).*(f1-f2).^2;
p = p./(max(p)+min(p)+eps);
end

% Eq. (20)
function [Ax,Ah,G,hG,CV]=Update_A_ordinary(Ax,CV,FE,G,hG,ub,lb,fun)
AubE = Ax>ub;
AlbE = Ax<lb;
Ax(:,AubE) = ub(AubE);
Ax(:,AlbE) = lb(AlbE);
Ah = fun (Ax);

if hG>Ah
    hG = Ah;
    G = Ax;
end
CV(FE) = hG;
end

% Eq. (20)
function [Bx,Bh,G,hG,CV]=Update_B_ordinary(Bx,CV,FE,G,hG,ub,lb,fun)
BubE = Bx>ub;
BlbE = Bx<lb;
Bx(:,BubE) = ub(BubE);
Bx(:,BlbE) = lb(BlbE);
Bh = fun (Bx);

if hG>Bh
    hG = Bh;
    G = Bx;
end
CV(FE) = hG;
end

% Eq. (19)
function [Ax,Ah,G,hG,CV]=Update_A_mod(Ax,CV,FE,G,hG,ub,lb,fun)
AubE = Ax>ub;
AlbE = Ax<lb;
Ax(:,AubE) = mod(Ax(:,AubE),ub(AubE)+eps)./(ub(AubE)+eps).*(ub(AubE)-lb(AubE)) + lb(AubE);
Ax(:,AlbE) = mod(Ax(:,AlbE),lb(AlbE)+eps)./(lb(AlbE)+eps).*(ub(AlbE)-lb(AlbE)) + lb(AlbE);
Ah = fun (Ax);

if hG>Ah
    hG = Ah;
    G = Ax;
end
CV(FE) = hG;
end

% Eq. (19)
function [Bx,Bh,G,hG,CV]=Update_B_mod(Bx,CV,FE,G,hG,ub,lb,fun)
BubE = Bx>ub;
BlbE = Bx<lb;
Bx(:,BubE) = mod(Bx(:,BubE),ub(BubE)+eps)./(ub(BubE)+eps).*(ub(BubE)-lb(BubE)) + lb(BubE);
Bx(:,BlbE) = mod(Bx(:,BlbE),lb(BlbE)+eps)./(lb(BlbE)+eps).*(ub(BlbE)-lb(BlbE)) + lb(BlbE);
Bh = fun (Bx);

if hG>Bh
    hG = Bh;
    G = Bx;
end
CV(FE) = hG;
end


% This function initialize the first population of search agents
function [x, new_lb, new_ub] = Initialization(N,dim,ub,lb)

    Boundary= size(ub,2); % numnber of boundaries
    new_lb = lb;
    new_ub = ub;
    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if Boundary==1
        x=rand(N,dim).*(ub-lb)+lb;
        new_lb = lb*ones(1,dim);
        new_ub = ub*ones(1,dim);
    end
    
    % If each variable has a different lb and ub
    if Boundary>1
        for i=1:dim
            ubi=ub(i);
            lbi=lb(i);
            x(:,i)=rand(N,1).*(ubi-lbi)+lbi;
        end
    end
end