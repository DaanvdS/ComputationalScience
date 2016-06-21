function [ out ] = d2Vdr2( r, sigm, eps )
    out=(624*sigm^12*eps*r.^(-14))-(168*eps*sigm^6*r.^(-8));
end