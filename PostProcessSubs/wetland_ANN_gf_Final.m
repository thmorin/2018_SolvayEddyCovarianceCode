% train an atificial Neural Network (ANN) to gap-filling
% missing or bad data of long term eddy-correlation measurements
%
% usage:
%  [dd_ann,holes,dd_ann2,tr_rec,weights1,weights2] = ANN_gf(t,p,arch,epochs,r2min)
% dd_ann   : gapfilled time series
% holes    :
% dd_ann2  : entirley reconstructed time series (useful for comparison with
%            original time series
% tr_rec   : Training record.
% weights1,weights2 : weights1 weights of first hidden layer
% arch    :  ANN architecture e.g. [12,5,1]
% t(n,1)  :  1-target (e.g., the flux series)
% p(n,m)  :  m-inputs (e.g., Rnet, taiir etc.)
% epochs  :  number of iterations
% r2min   :  minimum r2 between model and observations (if r2 is less than
%            r2min the training is repaeted from scratch and r2min is
%            reduced

function [dd_ann,holes,dd2_ann,r2max,tr_rec,weights1,weights2] = wetland_ANN_gf_Final(t,p,varargin)

%default options
opts=struct('arch',[12 5],...
    'epochs',100, ...
    'r2min',.75, ...
    'maxiter',25, ...
    'graph','n');
opts=parseArgs_gf(varargin,opts);
[trow, tcol]=size(t);
if trow>tcol
    t=t';
    p=p';
end
%paramenters and costants
use=find(isnan(p), 1);
if ~isempty(use); display('warning: p contains nans'); end

dd_ann=t;
dd2_ann=t;

use=find(~isnan(t));
gap=find(isnan(t));
%N=length(use);

%inizialize variables
r2_b=0;
iter=0;
r=1;
r2sum=0;
r2max=0;
offset=0.1;
while r2_b<opts.r2min*r %&& iter<opts.maxiter
    iter=iter+1;
    r=1-iter/100;
    
    if mod(iter,50)==0 && opts.r2min>r2sum/iter;
        opts.r2min=r2sum/iter + offset;
        offset=offset-0.02;
%         disp('Triggered!');
    end
    
    % Divide the dataset in two random blocks: one for training and one for testing
%    perm=randperm(N);
%    tran=use(perm(1:round(3*N/4)));
%    test=use(perm(round(N/4)+1:N));
    
    net = feedforwardnet(opts.arch);
    
    net.layers{1}.transferFcn = 'tansig';
    net.layers{2}.transferFcn = 'tansig';
    net.layers{3}.transferFcn = 'tansig';%;'purelin'
    net.trainFcn = 'trainlm';
    
    % Set up Division of Data for Training, Validation, Testing
    net.divideParam.trainRatio = 50/100;
    net.divideParam.valRatio = 25/100;
    net.divideParam.testRatio = 25/100;
    
    % normilize target and inputs between -1 and 1
    %[y1,PS] = mapminmax([p(:,tran);t(tran)]);
    [y1,PS] = mapminmax([p(:,use);t(use)]);
    p_use=y1(1:end-1,:);
    t_use=y1(end,:);
    
    % set taining parameters
    net.trainParam.show = NaN;
    net.trainParam.epochs = opts.epochs;
    net.trainParam.goal = 2.5e-04;
    net.trainParam.showWindow = false;
    
    %training the network
    [net,~,~] = adapt(net,p_use,t_use);
    [net,tr]=train(net,p_use,t_use,'useParallel','yes','showResources','no'); 
    
    norm_p=mapminmax.apply([p;nan(size(t))],PS); % Remap the whole driver matrix
    tann_=net(norm_p(1:end-1,:)); % Model t with the network found for the drivers matrix
    tann=mapminmax.reverse([nan(size(p)) ; tann_],PS); % Un-normal p and the output tann
    tann=tann(end,:); % Taks just the last row, that's were tann is
    
    plot(1:length(t),tann,1:length(t),t,'r')
    
    [R,Ptran]=corrcoef(t(use(tr.trainInd)),tann(use(tr.trainInd)));find(Ptran<0.05);
    r2_a=R(2,1)^2;r2sum=r2sum+r2_a;
    [R,Ptest]=corrcoef(t(use(tr.testInd)),tann(use(tr.testInd)));find(Ptest<0.05);
    r2_b=R(2,1)^2;
%     disp(['r2_b: ' num2str(r2_b)])

    % Ready for export out of this function if
    if r2_b>r2max
        dd_ann(gap)=tann(gap);
        dd2_ann=tann;
        r2max=r2_b;
        tr_rec=tr;
        weights1=net.IW;
        weights2=net.LW;
    end

end

%------------------------------------------------
% basic filtering
dd_ann(dd_ann<min(t) | dd_ann>max(t))=nan;
dd2_ann(dd2_ann<min(t) | dd2_ann>max(t))=nan;

if opts.graph=='y'
    figure(3)
    clf
    plot(t)
    hold on
    plot(dd2_ann,'.')
    % hold on
    % plot(decday,dd_ann,'o')
    legend('original data','ANN reconstructed series')
end

holes=nan;%dd_ann(isnan(t) & ~isnan(dd_ann));
