function [x,N_eval,N_iter,normg] = gaussnewton(phi,t,y,x0,tol,printout,plotout)
    x = x0;        
    r = @(z) (phi(z, t) - y);
    f = @(z) norm(r(z));
    N_eval = 0;
    N_iter = 1;
    epsilon = 0.00001;
    last_val = f(x);
    N_eval = N_eval + 1;
    
    fprintf("iter\tx\tmax(abs(r))\tnorm(grad)\tf evals\tlambda\n")

    while N_iter < 100


        J = jacobian(r, x);
        % for each direction we need 2 function evals
        N_eval = N_eval + 2 * numel(x0);
    
        d = -(J' * J + epsilon * eye(size(x0))) \ (J' * r(x));
        F = @(lambda) f(x + d * lambda);
        
        [lambda, evals] = linesearch(F, tol);
        N_eval = N_eval + evals;
        x = x + d * lambda;

        val = f(x);
        N_eval = N_eval + 1;
        if abs(val - last_val) < tol
            break;
        end
        last_val = val;

        if printout
            fprintf('%d\t%4.1f\t%4.4f\t\t%4.4f\t\t%d\t%4.4f\n', N_iter, x(1), max(abs(r(x))), norm(grad(f, x)), N_eval, lambda);
        end

        N_iter = N_iter + 1;
    end

    
    py = phi(x, t);

    if plotout
        plot(t, py);
        hold on;
        scatter(t, y, 'x');
        hold off;
    end

    normg = norm(grad(f, x));
end