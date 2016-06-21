function [ out ] = dVdr( r, sigm, eps )
    out=-(48*eps*sigm^12*r.^(-13))+(24*eps*sigm^6*r.^(-7));
end