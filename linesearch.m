function [lambda, N_evals] = linesearch(F, tol)
    %Start with armijo's rule
    alpha = 2;
    epsilon = 0.2;
    
    f0 = F(0);
    df0 = grad(F, 0);

    lambda = 0.1;

    N_evals = 4;

    while F(alpha * lambda) < f0 + epsilon * df0 * lambda * alpha
        lambda = alpha * lambda;
        N_evals = N_evals + 1;
    end
    
    N_evals = N_evals + 1;
    while F(lambda) > f0 + epsilon * df0 * lambda
        lambda = lambda / alpha;
        N_evals = N_evals + 1;
    end

%Then we switch to golden section search
    phi = (sqrt(5)-1)/2;
    a = 0;
    b = lambda * alpha;

    %simple loop to make sure F(b) is not inf
    while isinf(F(b))
        b = b / alpha;
        N_evals = N_evals + 1;
    end

    p1 = a + (1 - phi) * (b - a);
    f_1 = F(p1);

    p2 = a + phi * (b - a);
    f_2 = F(p2);
    N_evals = N_evals + 3;

    while abs(a-b) > tol
        if f_1 <= f_2
            b = p2;
            p2 = p1;
            f_2 = f_1;
            p1 = a + (1 - phi) * (b - a);
            f_1 = F(p1);
        else
            a = p1;
            p1 = p2;
            f_1 = f_2;
            p2 = a + phi * (b - a);
            f_2 = F(p2);
        end
        N_evals = N_evals + 1;
    end

    lambda = a;
    if isnan(F(lambda)) || F(lambda)>F(0)
        error('Bad job of the line search!')
    end
end