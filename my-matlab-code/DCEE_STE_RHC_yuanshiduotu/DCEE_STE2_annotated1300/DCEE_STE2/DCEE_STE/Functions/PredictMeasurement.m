function [Dsim] = PredictMeasurement(s,p,thresh)
%PredictMeasurement Summary of this function goes here
%   Detailed explanation goes here

Dsim = RadioactiveDispersionModel(s,p);
ersize = size(Dsim);
error = 0.05 * Dsim .* randn(ersize(1),ersize(2));
Dsim = Dsim + error;
Dsim(Dsim<thresh)=0;

if rand<0.3
    Dsim = 0;
end

end

