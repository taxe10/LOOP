function w = proj(u,v)
% This function projects vector v on vector u
w = (((dot(v,u)) /(dot(u,u))) * u);
end