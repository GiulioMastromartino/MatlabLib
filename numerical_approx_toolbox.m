function [res, err_est] = numerical_approx_toolbox(mode, func, interval, param)
%NUMERICAL_APPROX_TOOLBOX
%   mode: 'simpson', 'trapezi', 'interp_err'
%   param: N (intervalli) oppure h (passo per interp)

    a = interval(1); b = interval(2);
    err_est = NaN;
    
    switch mode
        case 'simpson'
            N = param; if mod(N,2)~=0, N=N+1; end
            h = (b-a)/N; x = linspace(a,b,N+1);
            y = func(x);
            w = ones(1,N+1); w(2:2:end-1)=4; w(3:2:end-2)=2;
            res = (h/3)*sum(w.*y);
            
        case 'trapezi'
            N = param; h = (b-a)/N; x = linspace(a,b,N+1);
            y = func(x);
            res = (h/2)*(y(1) + 2*sum(y(2:end-1)) + y(end));
            
        case 'interp_err'
            % Stima numerica max|f''|
            h_val = param;
            x_fine = linspace(a,b,1000); y_f = func(x_fine);
            d2 = gradient(gradient(y_f, x_fine(2)-x_fine(1)), x_fine(2)-x_fine(1));
            res = (h_val^2/8)*max(abs(d2));
    end
end
