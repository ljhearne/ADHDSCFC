function r = RosenthalR(z,N)
%Rosenthal, R. (1994). Parametric measures of effect size. In H. Cooper &
%L. V. Hedges (Eds.), The handbook of research synthesis. (pp. 231-244).
%New York: Russell Sage Foundation.

%effect size for ranksum t-tests
%r = Z/?N
r = z / (sqrt(N));
end

