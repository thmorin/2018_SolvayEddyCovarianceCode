function [Resp_obs,Resp_gf,Resp_Recon,Resp_Error,GPP_obs,GPP_gf,GPP_Recon,GPP_Error,NEE_obs,NEE_gf,NEE_reconstructed,NEE_Error,Resp_Full,GPP_Full,NEE_Full]=G_NN_OWRP(R,BadDataFlags,daynight,Dspk,Resp_obs,RespFillers,GPPFillers)
tic
%maybe make a toggle for resp in the winter? And leave some instructions on
%how to determine if it should be on or off

% Resp_gf=nan(size(daynight));
% Resp_Reconstructed=nan(size(daynight));
BestToKeep=ceil(.1*R);

Resp_gf=nan(length(daynight),1);
Resp_Recon=nan(length(daynight),1);
Resp_Error=nan(size(daynight));
Resp_gf_Full=nan(length(daynight),BestToKeep);
Resp_Full=nan(length(daynight),BestToKeep);

GPP_gf=nan(length(daynight),1);
GPP_Recon=nan(length(daynight),1);
GPP_gf_Full=nan(length(daynight),BestToKeep);
GPP_obs=nan(size(daynight));
GPP_Error=nan(size(daynight));
GPP_Full=nan(length(daynight),BestToKeep);

Resp_obs(max(BadDataFlags,[],2))=nan;
NEE_obs=Resp_obs;

Resp_obs(daynight)=nan;
Resp_obs=DeSpike(Resp_obs,Dspk.Interval,Dspk.STD,Dspk.Resp.min,Dspk.Resp.max,Dspk.trimp);

Resp_NN=nan(length(daynight),R);
Resp_Reconstructed_NN=nan(length(daynight),R);

%%%%%%%%%%%%%%%Respiration neural network
if sum(~isnan(Resp_obs))>20
    parfor i=1:R 
        [Resp_NN(:,i), ~, Resp_Reconstructed_NN(:,i),R2_Re(1,i)] = wetland_ANN_gf_Final(Resp_obs,RespFillers);
        [p,~,~,~,stats]=regress(Resp_obs,[Resp_Reconstructed_NN(:,i),ones(size(Resp_obs))]);
        disp(['Resp iteration:' num2str(i) ', r2 of ' num2str(stats(1)) ', slope of ' num2str(p(1))])
    end
    [~,I]=sort(R2_Re);

    Resp_Full=Resp_Reconstructed_NN(:,(I(end-BestToKeep+1:end)));
    Resp_gf_Full=Resp_NN(:,(I(end-BestToKeep+1:end)));
    Resp_gf=nanmean(Resp_gf_Full,2);
    Resp_Recon=nanmean(Resp_Reconstructed_NN(:,(I(end-BestToKeep+1:end))),2);
    
    GoodResp=~isnan(Resp_obs);
%    Resp_Error=sqrt(sum(Resp_obs(GoodResp)-Resp_Recon(GoodResp)).^2/sum(GoodResp))*ones(size(daynight));  
    Resp_Error=std(Resp_Reconstructed_NN,[],2);
    [p,~,~,~,stats]=regress(Resp_obs,[Resp_Recon,ones(size(Resp_obs))]);
    disp(['Overall NN Resp r2: ' num2str(stats(1)) ' slope of ' num2str(p(1))])
else
    disp('Warning: Less than 20 good respiration points. No attempt to gapfill was made.');
end



%%%%%%%%%%%GPP neural network
GPP_obs(daynight)=NEE_obs(daynight);
GPP_obs=GPP_obs-Resp_gf;
GPP_obs=DeSpike(GPP_obs,Dspk.Interval,Dspk.STD,Dspk.GPP.min,Dspk.GPP.max,Dspk.trimp);

GPP_NN=nan(sum(daynight),R);
GPP_Reconstructed_NN=nan(sum(daynight),R);

%GPP Loop for winter
if sum(~isnan(GPP_obs))>20
    parfor i=1:R 
        [GPP_NN(:,i), ~, GPP_Reconstructed_NN(:,i),R2_GPP(1,i)] = wetland_ANN_gf_Final(GPP_obs(daynight),GPPFillers(daynight,:));
        [p,~,~,~,stats]=regress(GPP_obs(daynight),[GPP_Reconstructed_NN(:,i),ones(size(GPP_obs(daynight)))]);
        disp(['GPP iteration:' num2str(i) ', r2 of ' num2str(stats(1)) ', slope of ' num2str(p(1))])
    end
    
    [~,I]=sort(R2_GPP);
    
    GPP_Full(daynight,:)=GPP_Reconstructed_NN(:,(I(end-BestToKeep+1:end)));
    GPP_Full(~daynight,:)=0;
    
    GPP_gf_Full(daynight,:)=GPP_NN(:,(I(end-BestToKeep+1:end)));
    GPP_gf_Full(~daynight,:)=0;
    GPP_gf=nanmean(GPP_gf_Full,2);
    GPP_Recon(daynight)=nanmean(GPP_Reconstructed_NN(:,(I(end-BestToKeep+1:end))),2);
    GPP_Recon(~daynight,:)=0;   
    GoodGPP=~isnan(GPP_obs);
%POP BREAK THE CODE SO I CANT RUN IT UNTIL I FIX STUFF
%	Tim: beef up the error stuff in here. Double check
%    GPP_Error=sqrt(sum(GPP_obs(GoodGPP)-GPP_Recon(GoodGPP)).^2/sum(GoodGPP))*ones(size(daynight));
    GPP_Error=std(GPP_Reconstructed_NN,[],2);
    [p,~,~,~,stats]=regress(GPP_obs(daynight),[GPP_Recon(daynight),ones(size(GPP_obs(daynight)))]);
    disp(['Overall GPP r2: ' num2str(stats(1)) ' slope of ' num2str(p(1))])
else
    disp('Warning: Less than 20 good GPP points. No attempt to gapfill was made');
end


%Combines neural networked respiration and GPP for NEE
NEE_obs=Resp_obs;
NEE_obs(isnan(NEE_obs))=Resp_gf(isnan(NEE_obs))+GPP_obs(isnan(NEE_obs));
NEE_gf=Resp_gf+GPP_gf;
NEE_Full=Resp_Full+GPP_Full;
NEE_reconstructed = Resp_Recon+GPP_Recon;
%NEE_Error=nansum([Resp_Error GPP_Error],2);
NEE_Error=std(NEE_Full,[],2);
toc

end
