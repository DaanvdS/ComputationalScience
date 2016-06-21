function [ out ] = dPdr( r, r_equi, ks )
    out=ks*(r-r_equi);
end