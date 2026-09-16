function [Qfunc] = illumination(sigma, F, tau1, tau2)
lim=tau2/tau1;
E=F*sigma;
b = 1 / tau1;
a = pi / tau2^2;
if tau2==Inf
        temp_fun = @(t) E * ones(size(t));
else
    if lim<.1
        temp_fun = @(t) E*b*exp(-b*t).*(t>=0);
        
    elseif lim>33
        delta=sqrt(log(1e3)/pi)*tau2;
        temp_fun = @(t) E/tau2*exp(-pi*(t-delta).^2/tau2^2);
    else
        delta=shift_convolution(tau1, tau2, 1e-3);
        temp_fun = @(t) E/(2*tau1) * exp(b^2 / (4 * a)) .* ...
            exp(-b * (t-delta)) .*(erfc(sqrt(a) * (-(t-delta) + b / (2 * a))));
           
    end
end
        % Función térmica pura (sin condición espacial)
        Qfunc = temp_fun;
end
