function [Qfunc] = illumination(sigma, F)

E=F*sigma;
temp_fun = @(t) E * ones(size(t));

% Heat source.
Qfunc = temp_fun;
end
