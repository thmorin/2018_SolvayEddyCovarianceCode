function [LE_obs,LE_gf,LE_Recon,LE_Error,H_obs,H_gf,H_Recon,H_Error,LE_Full,LE_gf_Full,H_Full,H_gf_Full]=G_NN_OWRP(R,BadDataFlags,daynight,LE_obs,H_obs,LEFillers,HFillers)
%Be sure that LE and H are despiked before feeding them into this function
tic
BestToKeep=ceil(.1*R);

LE_gf=nan(length(daynight),1);LE_Recon=nan(length(daynight),1);LE_Error=nan(size(daynight));LE_gf_Full=nan(length(daynight),BestToKeep);LE_Full=nan(length(daynight),BestToKeep);
H_gf=nan(length(daynight),1);H_Recon=nan(length(daynight),1);H_gf_Full=nan(length(daynight),BestToKeep);H_Error=nan(size(daynight));H_Full=nan(length(daynight),BestToKeep);

LE_obs(max(BadDataFlags,[],2))=nan;H_obs(max(BadDataFlags,[],2))=nan;
LE_NN=nan(length(daynight),R);LE_Reconstructed_NN=nan(length(daynight),R);

%%%%%%%%%%%%%%% LE neural network
if sum(~isnan(LE_obs))>20
    LE_NN_day=nan(sum(daynight),R);LE_Reconstructed_NN_day=nan(sum(daynight),R);
    LE_NN_night=nan(sum(~daynight),R);LE_Reconstructed_NN_night=nan(sum(~daynight),R);
    parfor i=1:R 
        [LE_NN_day(:,i), ~, LE_Reconstructed_NN_day(:,i),R2_LE_day(1,i)] = wetland_ANN_gf_Final(LE_obs(daynight),LEFillers(daynight,:),'r2min',0.9);
        [p,~,~,~,stats]=regress(LE_obs(daynight),[LE_Reconstructed_NN_day(:,i),ones(size(LE_obs(daynight)))]);
        disp(['LE day iteration:' num2str(i) ', r2 of ' num2str(stats(1)) ', slope of ' num2str(p(1))])
    end
    parfor i=1:R
        [LE_NN_night(:,i), ~, LE_Reconstructed_NN_night(:,i),R2_LE_night(1,i)] = wetland_ANN_gf_Final(LE_obs(~daynight),LEFillers(~daynight,:));
        [p,~,~,~,stats]=regress(LE_obs(~daynight),[LE_Reconstructed_NN_night(:,i),ones(size(LE_obs(~daynight)))]);
        disp(['LE night iteration:' num2str(i) ', r2 of ' num2str(stats(1)) ', slope of ' num2str(p(1))])
    end

    [~,I]=sort(R2_LE_day);
    LE_Full(daynight,:)=LE_Reconstructed_NN_day(:,(I(end-BestToKeep+1:end)));
    LE_gf_Full(daynight,:)=LE_NN_day(:,(I(end-BestToKeep+1:end)));

    [~,I]=sort(R2_LE_night);
    LE_Full(~daynight,:)=LE_Reconstructed_NN_night(:,(I(end-BestToKeep+1:end)));
    LE_gf_Full(~daynight,:)=LE_NN_night(:,(I(end-BestToKeep+1:end)));


    LE_gf=nanmean(LE_gf_Full,2);
    LE_Recon=nanmean(LE_Full,2);
    
    GoodLE=~isnan(LE_obs);
    LE_Error=sqrt(sum(LE_obs(GoodLE)-LE_Recon(GoodLE)).^2/sum(GoodLE))*ones(size(daynight));  
    [p,~,~,~,stats]=regress(LE_obs,[LE_Recon,ones(size(LE_obs))]);
    disp(['Overall NN LE r2: ' num2str(stats(1)) ' slope of ' num2str(p(1))])
else
    disp('Warning: Less than 20 good LE points. No attempt to gapfill was made.');
end

%%%%%%%%%%% H neural network
H_NN=nan(length(daynight),R);
H_Reconstructed_NN=nan(length(daynight),R);
if sum(~isnan(H_obs))>20
    H_NN_day=nan(sum(daynight),R);H_Reconstructed_NN_day=nan(sum(daynight),R);
    H_NN_night=nan(sum(~daynight),R);H_Reconstructed_NN_night=nan(sum(~daynight),R);
    parfor i=1:R 
        [H_NN_day(:,i), ~, H_Reconstructed_NN_day(:,i),R2_H_day(1,i)] = wetland_ANN_gf_Final(H_obs(daynight),HFillers(daynight,:),'r2min',0.9);
        [p,~,~,~,stats]=regress(H_obs(daynight),[H_Reconstructed_NN_day(:,i),ones(size(H_obs(daynight)))]);
        disp(['H day iteration:' num2str(i) ', r2 of ' num2str(stats(1)) ', slope of ' num2str(p(1))])
    end
    parfor i=1:R
        [H_NN_night(:,i), ~, H_Reconstructed_NN_night(:,i),R2_H_night(1,i)] = wetland_ANN_gf_Final(H_obs(~daynight),HFillers(~daynight,:));
        [p,~,~,~,stats]=regress(H_obs(~daynight),[H_Reconstructed_NN_night(:,i),ones(size(H_obs(~daynight)))]);
        disp(['H night iteration:' num2str(i) ', r2 of ' num2str(stats(1)) ', slope of ' num2str(p(1))])
    end
    
    [~,I]=sort(R2_H_day);   
    H_Full(daynight,:)=H_Reconstructed_NN_day(:,(I(end-BestToKeep+1:end)));
    H_gf_Full(daynight,:)=H_NN_day(:,(I(end-BestToKeep+1:end)));
    
    [~,I]=sort(R2_H_night);
    H_Full(~daynight,:)=H_Reconstructed_NN_night(:,(I(end-BestToKeep+1:end)));
    H_gf_Full(~daynight,:)=H_NN_night(:,(I(end-BestToKeep+1:end)));

    H_gf=nanmean(H_gf_Full,2);
    H_Recon=nanmean(H_Full,2);

    GoodH=~isnan(H_obs);
    H_Error=sqrt(sum(H_obs(GoodH)-H_Recon(GoodH)).^2/sum(GoodH))*ones(size(daynight));
    [p,~,~,~,stats]=regress(H_obs,[H_Recon,ones(size(H_obs))]);
    disp(['Overall H r2: ' num2str(stats(1)) ' slope of ' num2str(p(1))])
else
    disp('Warning: Less than 20 good H points. No attempt to gapfill was made');
end
toc
end
