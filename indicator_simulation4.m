function A=indicator_simulation4(cumprob_one)
    u_vec=rand([size(cumprob_one,1),1]);
    A=sum((cumprob_one-u_vec)<0,2)+1;
end