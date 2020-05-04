function finalpredict_Y = probDist(pred_Y,Ntest)

finalpredict_Y = zeros(Ntest,1);

for i = 1:Ntest
    Y = pred_Y(i,:);
    Y(isnan(Y))=[];
    a = unique(Y);
    if a == 0
        a = 1;
    end
    [x,y] = hist(Y,a);
    [~,indx] = max(x);
    finalpredict_Y(i,1) = y(indx);
end




end