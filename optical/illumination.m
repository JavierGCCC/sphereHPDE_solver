function [Qfunc] = illumination( F, tau1, tau2)

lim = tau2/tau1;   %Reduced timescale. 
E   = F;     %Absorbed power. 
b   = 1 / tau1;    %Characteristic factor.
a   = pi / tau2^2; %Characteristic factor.

if tau2==Inf %Check CW condition.
        temp_fun = @(t) E * ones(size(t)); %Constant input.
else %Check validity ranges.
    if lim<.1 %Short pulses.
        temp_fun = @(t) E*b*exp(-b*t).*(t>=0);
        
    elseif lim>33 %Long pulses.
        delta=sqrt(log(1e3)/pi)*tau2;
        temp_fun = @(t) E/tau2*exp(-pi*(t-delta).^2/tau2^2);
    else %Intermediate regimes.
        delta=shift_convolution(tau1, tau2, 1e-3); % Shift origin.
        temp_fun = @(t) E/(2*tau1) * exp(b^2 / (4 * a)) .* ...
            exp(-b * (t-delta)) .*...
            (erfc(sqrt(a) * (-(t-delta) + b / (2 * a))));
           
    end
end
        Qfunc = temp_fun; %Homogeneous distribution.
end
