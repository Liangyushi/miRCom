function [miRNAfun] = dmFunMiRNASim(meshSim,dm)
    %meshSim= importdata(f1,'\t'); 
    %DD = meshSim;
    dmsize=size(dm,2);
    cell_dm=cell(dmsize,1);
    for i=1:dmsize
        cell_dm{i}=find(dm(:,i)==1);
    end
    %nm=size(dm,2)
    %miRNAfun = zeros(541,541);
    miRNAfun = zeros(dmsize,dmsize);
    for i=1:dmsize
        miRNAfun(i,i)=1;
        Gi=cell_dm{i};
        len1=size(Gi,1);
        if len1~=0
            for j=1:dmsize
                if j ~= i
                    Gj=cell_dm{j};
                    len2=size(Gj,1);
                    if len2 ~=0
                        SGij=meshSim(Gi,Gj);
                        sumGi=sum(max(SGij,[],2)); %按行求
                        sumGj=sum(max(SGij,[],1)); %按列求
                        miRNAfun(i,j)= (sumGi+sumGj)/(len1+len2);
                        miRNAfun(j,i)=miRNAfun(i,j);
                    end
                end
            end
        end
    end
                
end

