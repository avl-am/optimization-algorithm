%[2017]-"A new meta-heuristic butterfly-inspired algorithm"Artificial
%Butterfly Optimization - ABO

function [bestf,bestx,cg_curve] = ABO(N, max_Iter, lb,ub ,dim, fun)
% Parameters

step_e = 0.05;   % control number of sunspot
ratio  = 0.2;    % control step
type   = 1;      % type 1 or 2
% Initial
X   = zeros(N,dim);
for i = 1:N
    % 	for d = 1:dim
    X(i,:) = lb + (ub - lb).* rand();
    %   end
end
% Fitness
fit  = zeros(1,N);
fitG = inf;
for i = 1:N
    fit(i) = fun(X(i,:));
    % Global update
    if fit(i) < fitG
        fitG = fit(i);
        Xgb  = X(i,:);
    end
end
% Pre
Xnew = zeros(N,dim);
curve = zeros(1,max_Iter);
curve(1) = fitG;
t = 2;
% Iteration
while t <= max_Iter
    % Sort butterfly
    [fit, idx] = sort(fit,'ascend');
    X          = X(idx,:);
    % Proportion of sunspot butterfly decreasing from 0.9 to ratio
    num_sun = round(N * (0.9 - (0.9 - ratio) * (t / max_Iter)));
    % Define a, linearly decrease from 2 to 0
    a       = 2 - 2 * (t / max_Iter);
    % Step update (5)
    step    = 1 - (1 - step_e) * (t / max_Iter);
    % {1} Some butterflies with better fitness: Sunspot butterfly
    for i = 1:num_sun
        % Random select a butterfly k, but not equal to i
        R = randperm(N); R(R == i) = [];
        k = R(1);
        % [Version 1]
        if type == 1
            % Randomly select a dimension
            J  = randi([1,dim]);
            % Random number in [-1,1]
            r1 = -1 + 2 * rand();
            % Position update (1)
            Xnew(i,:) = X(i,:);
            Xnew(i,J) = X(i,J) + (X(i,J) - X(k,J)) * r1;
            % [Version 2]
        elseif type == 2
            % Distance
            dist = norm(X(k,:) - X(i,:));
            r2   = rand();
            for d = 1:dim
                % Position update (2)
                Xnew(i,d) = X(i,d) + ((X(k,d) - X(i,d)) / dist) * ...
                    (ub - lb) * step * r2;
            end
        end
        % Boundary
        XB = Xnew(i,:);
        %     XB(XB > ub) = ub;
        %     XB(XB < lb) = lb;
        
        Flag4ub=XB>ub;
        Flag4lb=XB<lb;
        XB=(XB.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        Xnew(i,:) = XB;
    end
    % Fitness
    for i = 1:num_sun
        % Fitness
        Fnew = fun(X(i,:));
        % Greedy selection
        if Fnew < fit(i)
            fit(i) = Fnew;
            X(i,:) = Xnew(i,:);
        end
        % Global update
        if fit(i) < fitG
            fitG = fit(i);
            Xgb  = X(i,:);
        end
    end
    
    % {2} Some butterflies: Canopy butterfly
    for i = num_sun + 1 : N
        % Random select a sunspot butterfly
        k = randi([1,num_sun]);
        % [Version 1]
        if type == 1
            % Randomly select a dimension
            J  = randi([1,dim]);
            % Random number in [-1,1]
            r1 = -1 + 2 * rand();
            % Position update (1)
            Xnew(i,:) = X(i,:);
            Xnew(i,J) = X(i,J) + (X(i,J) - X(k,J)) * r1;
            % [Version 2]
        elseif type == 2
            % Distance
            dist = norm(X(k,:) - X(i,:));
            r2   = rand();
            for d = 1:dim
                % Position update (2)
                Xnew(i,d) = X(i,d) + ((X(k,d) - X(i,d)) / dist) * ...
                    (ub - lb) * step * r2;
            end
        end
        % Boundary
        XB = Xnew(i,:);
        
        %         XB(XB > ub) = ub; XB(XB < lb) = lb;
        Flag4ub=XB>ub;
        Flag4lb=XB<lb;
        XB=(XB.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        
        Xnew(i,:) = XB;
    end
    % Fitness
    for i = num_sun + 1 : N
        % Fitness
        Fnew = fun(X(i,:));
        % Greedy selection
        if Fnew < fit(i)
            fit(i) = Fnew;
            X(i,:) = Xnew(i,:);
        else
            % Random select a butterfly
            k  = randi([1,N]);
            % Fly to new location
            r3 = rand();
            r4 = rand();
            for d = 1:dim
                % Compute D (4)
                Dx     = abs(2 * r3 * X(k,d) - X(i,d));
                % Position update (3)
                X(i,d) = X(k,d) - 2 * a * r4 - a * Dx;
            end
            % Boundary
            XB = X(i,:);
            
            %             XB(XB > ub) = ub; XB(XB < lb) = lb;
            Flag4ub=XB>ub;
            Flag4lb=XB<lb;
            XB=(XB.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            
            
            X(i,:) = XB;
            % Fitness
            fit(i) = fun(X(i,:));
        end
        % Global update
        if fit(i) < fitG
            fitG = fit(i);
            Xgb  = X(i,:);
        end
    end
    curve(t) = fitG;
    t = t + 1;
end
bestf = curve(end);
cg_curve = curve;
bestx = Xgb;
end
