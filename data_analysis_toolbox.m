function [out1, out2, out3] = data_analysis_toolbox(mode, arg1, arg2, arg3)
%DATA_ANALYSIS_TOOLBOX
%   mode: 'lsq' (Minimi Quadrati), 'roots' (Analisi Newton), 'svd' (Approx Rango)
    
    switch mode
        case 'lsq'
            % arg1=x, arg2=y, arg3=degree
            coeffs = polyfit(arg1, arg2, arg3);
            y_fit = polyval(coeffs, arg1);
            resid_sq = sum((arg2 - y_fit).^2);
            out1 = coeffs; out2 = resid_sq; out3 = [];
            
        case 'roots'
            % arg1 = errors array
            errs = arg1;
            if errs(1)<errs(end), errs=abs(errs(1:end-1)-errs(end)); end % Convert to errors
            
            % Stima p
            ratios = log(errs(3:end)./errs(2:end-1)) ./ log(errs(2:end-1)./errs(1:end-2));
            p_est = ratios(end);
            C_est = errs(end)/(errs(end-1)^p_est);
            m_est = 1/(1-C_est); if abs(p_est-1)>0.5, m_est=1; end
            
            out1 = p_est; out2 = C_est; out3 = m_est;
            
        case 'svd'
            % arg1=A, arg2=k
            [U,S,V] = svd(arg1);
            diag_S = diag(S);
            
            Uk = U(:,1:arg2); Sk = S(1:arg2, 1:arg2); Vk = V(:,1:arg2);
            Ak = Uk*Sk*Vk';
            err_rel = 0; if arg2<length(diag_S), err_rel = diag_S(arg2+1)/diag_S(1); end
            
            out1 = Ak; out2 = err_rel; out3 = diag_S;
    end
end
