function dy =hydrogenODE_P(t,y,m_eff,ks,r_equi)
dy = [y(2) ; -ks*(y(1)-r_equi)/m_eff];
end