function resultSet=SubCodeNormRSDMS(dataPack,globKMat,nodeLabelList,constrainedDoFVec,elaMatrix,numSet,degree,tolSet)

elemYield=cellfun(@(x)x.Yield,dataPack,'UniformOutput',0);
elaStress=getElemGaussStress(dataPack);  % Get Gauss-point elastic stress
degreeSet=deg2rad(degree);
weightSet=[cos(degreeSet),sin(degreeSet),zeros(numel(degreeSet),1)];

resultSet=cell(numel(degreeSet),1);      % Preallocate result container

%% 定义路径方法1：
[loadingCoe,GaussPoint,GaussWeight,numSet]=subLoadingCoe(numSet);

%% 定义路径方法2
% [starredPoint,GaussPoint,GaussWeight]=getStarredPoint(numSet);
% [path1,path2,numSet]=subLoadingPathFunc(numSet);
% loadingCoe=[path1(starredPoint),path2(starredPoint)];
% loadingCoe=transpose(loadingCoe);

%%
timePoint=GaussPoint;  % time points (Gauss–Legendre points)

%% Visualize the selected loading path
% subLoadingPathVisualize(numSet,timePoint);

%%
elaFactor=1;
eleStressFull=elaFactor.*[elaStress,zeros(size(elaStress,1),3-size(elaStress,2))];

for iter=1:numel(degreeSet)
    
    tt=tic;  % Start timer
    weight=weightSet(iter,1:numSet.numLoadCase);
    vertWeight=diag(weight)*loadingCoe(1:numSet.numLoadCase,:);
    vertElaStress=eleStressFull(:,1:numSet.numLoadCase)*vertWeight;

    %% 初始化安定极限载荷乘子: 计算各高斯点的等效弹性应力,计算弹性极限;
    elemElaStress=mat2cell(vertElaStress,numSet.numStrComp*numSet.NGK*ones(numSet.numElem,1),numSet.numGL);
    elemElaStress=cellfun(@(x)mat2cell(x,numSet.numStrComp*ones(numSet.NGK,1),ones(1,numSet.numGL)),elemElaStress,'UniformOutput',0);
    maxElemEffStress=cellfun(@(x)max(max(cellfun(@subCalEffStress,x))),elemElaStress);
    elasticLimit=min(cell2mat(elemYield)./maxElemEffStress);
    alpha=numSet.alphaScaling*elasticLimit;  

    %% 初始化残余应力场 rho(t)
    elemRsiStress=mat2cell(sparse(numSet.NGK*numSet.numStrComp,1),numSet.numStrComp*ones(numSet.NGK,1),1);
    elemRsiStress=repmat({repmat(elemRsiStress,1,numSet.numGL)},numSet.numElem,1);
    
    %% 初始化傅里叶系数 a0 和 ak：rho(t) = 1/2 * a0 + sum_k [ak * cos(2kπt) + bk * sin(2kπt)]
    a0Coe=repmat({zeros(numSet.NGK,numSet.numStrComp)},numSet.numElem,1);
    akCoe=cellfun(@(x)repmat({x},numSet.numTerms,1),a0Coe,'UniformOutput',0);
    % bkCoe implicitly handled later (not initialized here)
    
    %% 初始化 cos(2kπt) 和 sin(2kπt)
    numTermsCell=num2cell(transpose(1:numSet.numTerms));
    cosComp=cellfun(@(t)cellfun(@(k)cospi(2*k*t),numTermsCell),num2cell(timePoint),'UniformOutput',0);
    sinComp=cellfun(@(t)cellfun(@(k)sinpi(2*k*t),numTermsCell),num2cell(timePoint),'UniformOutput',0);
    
    % Logging variables
    iterSet=[];
    alphaSet=[];
    deltaLambda1Set={};
    coePhiAlpha=[];
    coePhi=0;
    outerLoop=0;
    
    %% Outer iteration loop
    while true
        
        t=tic;
        innerLoop=0;
        outerLoop=outerLoop+1;
        alphaSet(outerLoop,1)=alpha; %#ok<AGROW>
        deltaLambda1=zeros(numSet.NGK*numSet.numElem,numel(timePoint));
        deltaLambda1=mat2cell(deltaLambda1,numSet.NGK*ones(numSet.numElem,1),numel(timePoint));

        %% Inner iteration loop (gradient descent-like)
        while true

            innerLoop=innerLoop+1;
            coePhi_old=coePhi;
            stress=cellfun(@(x,y)cellfun(@(e,r)alphaSet(end).*e+r,x,y,'UniformOutput',0),elemElaStress,elemRsiStress,'UniformOutput',0);
            %%%由总应力截断的超出屈服面部分作为塑性应力
            [plaStress, plaCoe]=cellfun(@(x,y)subCalElemPlaStress(x,y),stress,elemYield,'UniformOutput',0);
            %%%由总应力通过应力投影求解的垂直屈服面部分作为塑性应力
            % [deltaLambda1,plaStress]=cellfun(@(x,l,t,p)subCalLambda1Func(x,l,t,p),dataPack,deltaLambda1,stress,plaStress,'UniformOutput',0);
            % deltaLambda1=cellfun(@(x,l,t,p,c)subCalLambda1Func(x,l,t,p,c),dataPack,deltaLambda1,stress,plaStress,plaCoe,'UniformOutput',0);
            deltaLambda1=cellfun(@(x,l,t,p)subCalLambda1Func(x,l,t,p),dataPack,deltaLambda1,stress,plaStress,'UniformOutput',0);
            resiNodeForce=cellfun(@(x,p)subCalResiNodeForce_RSDMS(x,p),dataPack,plaStress,'UniformOutput',0);
            globResiNodeForce=subSortNodeForce2Model_RSDMS(dataPack,resiNodeForce,numel(nodeLabelList));
            globResiDisp=subCalGlobResiDisp(globResiNodeForce,globKMat,constrainedDoFVec);
            resiDisp=cellfun(@(x)subSortGlobResiDisp2Elem(x,globResiDisp),dataPack,'UniformOutput',0);
            totalResiStressRate=cellfun(@(x,y)subCalTotalResiStressRate(x,y,elaMatrix),dataPack,resiDisp,'UniformOutput',0);
            resiStressRate=cellfun(@(x,y)subResiDerivativeRate(x,y),totalResiStressRate,plaStress,'UniformOutput',0);
            [a0Coe,akCoe,bkCoe]=subUpdateBasisCoe(a0Coe,akCoe,timePoint,GaussWeight,resiStressRate);            
            [coePhi,coePhiTol]=UpdateCoePhionPlaStress(coePhi_old,plaStress,GaussWeight);
            elemRsiStress=cellfun(@(a0,a,b)subUpdateResiStress_RSDMS(a0,a,b,cosComp,sinComp),a0Coe,akCoe,bkCoe,'UniformOutput',0);
            if (coePhiTol<=tolSet.innerTol)||(innerLoop>=10)
                break;
            end

        end
        deltaLambda1=subNormUpdate(deltaLambda1,stress,elemElaStress);
        deltaLambda1Set{outerLoop,1}=cell2mat(deltaLambda1);
        coePhiAlpha(outerLoop,1)=coePhi;
        iterSet(outerLoop,1)=innerLoop;
        alpha=subUpdateAlpha_RSDMS(alphaSet,coePhiAlpha,outerLoop);
        %% Print iteration info
        % iterInfoOutput(outerLoop, innerLoop, alphaSet, coePhiAlpha, t);
        
        if coePhi <= tolSet.outerTol % Outer convergence check
            break;
        end

    end
    %% Post-processing and result storage
    tt=toc(tt);
    result=subPEMPostProcess(stress,elemRsiStress,elemYield);
    result.numSet=numSet;
    result.alpha=alpha;
    result.alphaSet=alphaSet;
    result.coePhiAlpha=coePhiAlpha;
    result.timePoint=timePoint;
    result.degree=degreeSet(iter);  % Angle in degrees
    result.vertWeight=vertWeight;
    result.vertElaStress=vertElaStress;
    result.deltaLambda1Set=deltaLambda1Set;

    Lambda1=sum(cat(3,deltaLambda1Set{end-3:end}),3);
    % Lambda1=sum(cat(3,deltaLambda1Set{end-3:end}),3);
    Lambda1=mat2cell(Lambda1,numSet.NGK*ones(numSet.numElem,1),size(Lambda1,2));
    Lambda1=subNormUpdate(Lambda1,stress,elemElaStress);

    result.Lambda1=Lambda1;
    result.work=tt;
    resultSet{iter}=result;

    %% Optional: plot shakedown limit curve
    if numel(degreeSet)>1
        plotShakedownDomain(resultSet,elemYield{1});  % Optional: plot shakedown limit curve
    end

end

end
