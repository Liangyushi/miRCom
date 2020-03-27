function [auc,fpr_tpr,aupr,r_p,finalmetrice] = Metrics(testclass,scores)
    [m,~] = size(testclass);
    len=size(testclass,1);
    tpr = zeros(1,len);
    fpr = zeros(1,len);
    recall = zeros(1,len);
    precision = zeros(1,len);
    f1 = zeros(1,len);
    num = 0;
    for x = 1:len
        num=num+1;
        tp = 0;
        fn = 0;
        fp = 0;
        tn = 0;
        for y = 1:m
            if scores(y,1)>= scores(x,1) && testclass(y,1)== 1
                tp=tp+1;
            elseif scores(y,1)>=scores(x,1) && testclass(y,1)==0
                fp=fp+1;
            elseif scores(y,1)<scores(x,1) && testclass(y,1)==0
                tn=tn+1;
            elseif scores(y,1)<scores(x,1) && testclass(y,1)==1
                fn=fn+1;
            end
        end
        tpr(1,num) = tp/(tp+fn);
        fpr(1,num) = fp/(fp+tn);
        recall(1,num)=tpr(1,num); 
        precision(1,num)=tp/(tp+fp); 
        f1(1,num)=2*((recall(1,num)*precision(1,num))/(recall(1,num)+precision(1,num)));
    end

    fpr_tpr = [fpr;tpr];
    r_p=[recall;precision];
    rfinal=sum(recall)/len;
    pfinal=sum(precision)/len;
    F1final=sum(f1)/len;

    
auc = -trapz(fpr,tpr);
aupr= -trapz(recall,precision);
finalmetrice=[auc aupr rfinal pfinal F1final];

end
