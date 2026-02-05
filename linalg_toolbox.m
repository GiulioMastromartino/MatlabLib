function [res, info] = linalg_toolbox(mode, A, opts)
%LINALG_TOOLBOX
%   mode: 'is_def_pos', 'power_iter', 'inv_power_iter'

    switch mode
        case 'is_def_pos'
            try, chol(A); res=true; catch, res=false; end
            info.eigs = eig(A);
            
        case 'power_iter'
            x = opts.x0(:); k = opts.k;
            for i=1:k
                y = A*x;
                [~,idx] = max(abs(y));
                res = y(idx); % Lambda approx
                x = y/norm(y, inf);
            end
            info.vec = x;
            
        case 'inv_power_iter'
            x = opts.x0(:); k = opts.k; 
            s = 0; if isfield(opts,'shift'), s=opts.shift; end
            [L,U,P] = lu(A - s*eye(size(A)));
            for i=1:k
                y = U\(L\(P*x));
                x = y/norm(y);
                res = (x'*A*x)/(x'*x); % Rayleigh
            end
            info.vec = x;
    end
end
