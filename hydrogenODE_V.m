function dy =hydrogenODE_V(t,y,m_eff,eps,sigm)
dy = [y(2) ; (48*eps*sigm^12*y(1)^(-13))/m_eff-(24*eps*sigm^6*y(1)^(-7))/m_eff];
end