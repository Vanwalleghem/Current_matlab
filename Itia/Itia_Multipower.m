Stimuli=zeros(6,size(ZS,2));
GCaMP6=[0,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
idxStart=58;
counter=0;
for i=1:18
    if mod(i,3)==1
        counter=counter+1;
    end
    Stimuli(counter,(idxStart+(i-1)*40):(idxStart+(i-1)*40)+size(GCaMP6,1)-1)=GCaMP6;
end

ModelMultipower=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=fitlm(Stimuli',ZS(i,:),'interactions');
    ModelMultipower(i).coef=mdl.Coefficients;
    %ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    %ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelMultipower(i).rsquared=mdl.Rsquared.Adjusted;
end

coefficients={};
for idx=1:length(ModelMultipower)
    coef=[ModelMultipower(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d+)','tokens');
    if ~isempty(temp)
        %temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=2:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp{coef_idx}{1}{1})}=coef.Estimate(coef_idx);
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients);
coefficients(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients=cell2mat(coefficients);