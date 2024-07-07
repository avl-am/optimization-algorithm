function [gBestScore, gBest, GlobalBestCost, AverageCost, StandardDCost] = CPSOGSA(n, iteration, low, up, dim, fobj)
    phi1 = 2.05;
    phi2 = 2.05;
    phi = phi1 + phi2;
    chi = 2 / (phi - 2 + sqrt(phi^2 - 4 * phi));
    w = chi;
    wdamp = 1;
    C1 = chi * phi1;
    C2 = chi * phi2;

    current_fitness = zeros(n, 1);
    gBest = zeros(1, dim);
    gBestScore = inf;
    VelMax = 0.1 * (up - low);
    VelMin = -up;

    for i = 1:n
        pBestScore(i) = inf;
    end
    pBest = zeros(n, dim);

    G0 = 1;
    current_position = rand(n, dim) .* (up - low) + low;
    velocity = 0.3 * randn(n, dim);
    acceleration = zeros(n, dim);
    mass = zeros(1, n);
    force = zeros(n, dim);

    iter = 0;
    GlobalBestCost = zeros(1, iteration);
    AverageCost = zeros(1, iteration);
    StandardDCost = zeros(1, iteration);

    while (iter < iteration)
        G = G0 * exp(-23 * iter / iteration);
        iter = iter + 1;

        for i = 1:n
            Tp = current_position(i, :) > up;
            Tm = current_position(i, :) < low;
            current_position(i, :) = (current_position(i, :) .* ~(Tp + Tm)) + up .* Tp + low .* Tm;

            fitness = fobj(current_position(i, :));
            current_fitness(i) = fitness;

            if (pBestScore(i) > fitness)
                pBestScore(i) = fitness;
                pBest(i, :) = current_position(i, :);
            end

            if (gBestScore > fitness)
                gBestScore = fitness;
                gBest = current_position(i, :);
            end
        end

        best = min(current_fitness);
        AVE = mean(current_fitness);
        STD = std(current_fitness);
        worst = max(current_fitness);

        GlobalBestCost(iter) = gBestScore;
        AverageCost(iter) = AVE;
        StandardDCost(iter) = STD;

        for pp = 1:n
            if current_fitness(pp) == best
                break;
            end
        end

        bestIndex = pp;

        for pp = 1:dim
            best_fit_position(iter, 1) = best;
            best_fit_position(iter, pp + 1) = current_position(bestIndex, pp);
        end

        mass = (current_fitness - 0.99 * worst) / (best - worst);
        mass = mass * 5 / sum(mass);

        force = zeros(n, dim);

        for i = 1:n
            for j = 1:dim
                for k = 1:n
                    if (current_position(k, j) ~= current_position(i, j))
                        force(i, j) = force(i, j) + rand() * G * mass(k) * mass(i) * (current_position(k, j) - current_position(i, j)) / abs(current_position(k, j) - current_position(i, j));
                    end
                end
            end
        end

        acceleration = zeros(n, dim);

        for i = 1:n
            for j = 1:dim
                if (mass(i) ~= 0)
                    acceleration(i, j) = force(i, j) / mass(i);
                end
            end
        end

        for i = 1:n
            for j = 1:dim
                velocity(i, j) = w * rand() * velocity(i, j) + C1 * rand() * acceleration(i, j) + C2 * rand() * (gBest(j) - current_position(i, j));
            end
        end

        current_position = current_position + velocity;
    end
end
