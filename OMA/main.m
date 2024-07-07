%% Optical Microscope Algorithm (OMA)
%  Author and programmer:
%         Professor        Min-Yuan Cheng
%         Ph.D. Student    Moh Nur Sholeh
%  Written by Moh Nur Sholeh
%  Computer Integrated Construction (CIC) Lab
%  National Taiwan University of Science and Technology, Taipei, Taiwan
%  Paper : Cheng, M. Y., & Sholeh, M. N. (2023). Optical microscope algorithm: a new metaheuristic inspired by microscope magnification
%          for solving engineering optimization problems. Knowledge-Based Systems.
%  DOI   : https://doi.org/10.1016/j.knosys.2023.110939 
%---------------------------------------------------------------------------------------------------------------------------%

function main()
    clear all
    clc
    tic;

    %% Select function
    funnum = 1;         % Function number
    [lb, ub, dim] = boundcondition(funnum);

    %% Set the parameters
    NP = 50;            % Number of target objects
    NI = 1000;          % Number of iterations
    para = [NI, NP];    % Parameter matrix

    %% Run OMA
    tic;
    [bestSolution, bestMagnification, numEval, bestMagnificationHistory] = OMA(@fobj, funnum, lb, ub, dim, para);
    time = toc;

    %% Display optimal results
    disp('OPTICAL MICROSCOPE ALGORITHM (OMA) RESULTS');
    disp('------------------------------------------');
    disp(['The best solution found by the OMA is: ' num2str(bestSolution)]);
    disp(['The highest magnification value achieved for the objective function is: ' num2str(bestMagnification)]);

    toc;

    %% Save optimal results
    save('result.mat', 'time', 'bestSolution', 'bestMagnification', 'numEval', 'bestMagnificationHistory');
end

function [lb, ub, dim] = boundcondition(funnum)
    % Define the bounds and dimension based on the function number
    switch funnum
        case 1
            lb = [-5 -5];
            ub = [5 5];
            dim = 2;
        case 2
            % Define bounds and dimension for function 2
        % Add more cases for other functions
    end
end

function y = fobj(x, funnum)
    % Define your objective functions here
    switch funnum
        case 1
            y = x(1)^2 + x(2)^2;  % Example objective function
        case 2
            % Define objective function 2
        % Add more cases for other objective functions
    end
end
