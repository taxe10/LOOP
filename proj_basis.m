function w = proj_basis(e,v)
% This function projects vector v on vector ui (e = u1,u2,...)
w = 0;
for i=1:size(e,1)
    w = w + (((dot(v,e(i,:))) /(dot(e(i,:),e(i,:)))) * e(i,:));
end
end