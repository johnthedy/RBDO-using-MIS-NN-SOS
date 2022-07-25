function [Pf,FE,firstbeta,DATASET,DATATESTING]=fobj2(eco,dring1,dring2,nagent1_1,nagent1_2,covMCS,nagent2,stopval,mu,sigma,nRV,dist,Rub,Rlb,AIstatus,net)
sphere=ceil(nRV*1.5);
counternum=5;
linearfold1=4;
linearfold2=4;
DATASET=[];
DATATESTING=[];

%% Start
candidatepoint=zeros(1,nRV);

ucenter=candidatepoint(1,:);
[sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,ucenter);
centerFE=G(eco,sample);FE=1;
if centerFE<0
    category=0;
else
    category=1;
end
DATASET=vertcat(DATASET,[sample category]);

if centerFE<0
    Pf=0.5;
    firstbeta=0;
else
    %beta1 is outter beta and beta2 is inner beta
    beta1=Rub;beta2=Rub-0.5;bestbeta=beta1;
    beta12=[];recdummy1_1=[];recdummy1_2=[];recfdummy1_2=[];totalbeta=[];stopcrit=0;tempFE=0;
    while beta1>0
        dummy0=0;
        p1=chi2cdf(beta1^2,nRV);
        p2=chi2cdf(beta2^2,nRV);
        beta12=vertcat(beta12,[beta1 beta2]);
        temp=sqrt(sum(recdummy1_2.^2,2));
        samplelimit=(1-(1-normcdf(min(beta1,3),0,1)))/((0.05^2)*(1-normcdf(min(beta1,3),0,1)))*(p1-p2)-length(temp((temp<beta1 & temp>beta2)));
        clear temp
        
        counter=counternum;
        while dummy0<nagent1_1
            if (tempFE>0.3*samplelimit && dummy0==0) || (tempFE>0.6*samplelimit && dummy0<2) || (tempFE>0.9*samplelimit && dummy0<3)
                break
            end
            x=normrnd(0,1,1,nRV);
            temp=unifrnd(p2,p1);
            x=x.*(sqrt(chi2inv(temp,nRV))/sqrt(sum(x.^2)));
            [sample,usample,~,~]=summonsample(1,mu,sigma,nRV,dist,x);DATATESTING=vertcat(DATATESTING,sample);
            if AIstatus==1 && counter<0
                label=round(net(sample'))';
                if label==0
                    dummy1=G(eco,sample);
                    if dummy1<0
                        category=0;
                    else
                        category=1;
                    end
                    recdummy1_1=vertcat(recdummy1_1,dummy1);recdummy1_2=vertcat(recdummy1_2,usample);
                    DATASET=vertcat(DATASET,[sample category]);
                    FE=FE+1;
                else
                    dummy1=1;
                end
            else
                dummy1=G(eco,sample);
                if dummy1<0
                    category=0;
                else
                    category=1;
                end
                DATASET=vertcat(DATASET,[sample category]);
                recdummy1_1=vertcat(recdummy1_1,dummy1);recdummy1_2=vertcat(recdummy1_2,usample);
                FE=FE+1;
                counter=counter-1;
            end
            tempFE=tempFE+1;
            if dummy1<0
                dummy0=dummy0+1;
                dummy1_1(dummy0,:)=dummy1;dummy1_2(dummy0,:)=usample;
            end
            if tempFE>samplelimit
                break
            end
        end
        nagent1_1=nagent1_2;
        clear dummy1 counter
        if dummy0==0
            break
        end
        for h=1:dummy0
            segmentc(1,:)=0.5;increment=0.5;
            for i=1:linearfold1
                increment=increment/2;
                dummy2(i,:)=(dummy1_2(h,:)-ucenter).*segmentc(i,:)+ucenter;
                dummy3(i,:)=sqrt(sum((dummy2(i,:)-ucenter).^2));
                [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,dummy2(i,:));DATATESTING=vertcat(DATATESTING,sample);
                dummy5(i,:)=G(eco,sample);
                if dummy5(i,:)<0
                    category=0;
                else
                    category=1;
                end
                DATASET=vertcat(DATASET,[sample category]);
                FE=FE+1;
                if dummy5(i,:)<0 && i<linearfold1
                    segmentc(i+1,:)=segmentc(i,:)-increment;
                else
                    segmentc(i+1,:)=segmentc(i,:)+increment;
                end
            end
            
            segmentc(end,:)=[];[~,idx]=sort(segmentc);
            segmentc=segmentc(idx);
            dummy3=dummy3(idx);
            dummy5=dummy5(idx);
            
            %Quadratic FIT
            dummy1=vertcat(centerFE,dummy5,dummy1_1(h));
            dummy2=vertcat(0,dummy3,beta1);
            dummy3=vertcat(0,segmentc,1);
            clear dummy4 dummy5 idx
            for i=1:length(dummy1)
                if dummy1(i)<0
                    break
                end
            end
            p=polyfit(dummy2(i-1:i),dummy1(i-1:i),1);
            storedbeta(h)=-p(2)/p(1);
            p=polyfit(dummy3(i-1:i),dummy1(i-1:i),1);
            dummy4=-p(2)/p(1);
            recfdummy1_2=vertcat(recfdummy1_2,(dummy1_2(h,:)-ucenter).*dummy4+ucenter);
            clear dummy1 dummy2 dummy3 dummy4 segmentc p
        end
        clear h i j
        
        if min(storedbeta)<beta2
            bestbeta=min(storedbeta);
            
            beta1=min(storedbeta);
            beta2=min(storedbeta)-dring1;
            stopcrit=0;
            tempFE=0;
        else
            stopcrit=stopcrit+1;
        end
        if stopcrit==stopval || beta1<=Rlb
            break
        end
    end
    totalbeta=vertcat(totalbeta,bestbeta);
    
    if totalbeta==Rub
        Pf=0;
        firstbeta=inf;
    elseif totalbeta<=Rlb
        Pf=0.5;
        firstbeta=0;
    else
        %% Candidate Point
        [candidatepoint1,~]=cp(totalbeta(1),recdummy1_1,recdummy1_2,nRV);
        clear beta1 beta2 dummy0 dummy1_1 dummy1_2 p1 p2 storedbeta tempFE temp lim x sample usample centerFE
        
        recdummy1_3=sqrt(sum(recdummy1_2.^2,2));
        betasegment=unique(beta12);
        idx=find(recdummy1_3<betasegment(1));recdummy1_1(idx,:)=[];recdummy1_2(idx,:)=[];
        clear beta12 idx
        
        %% Generate Sample
        Nsample=(1-(1-normcdf(totalbeta(1),0,1)))/((covMCS^2)*(1-normcdf(totalbeta(1),0,1)));
        if Nsample>(1-(1-normcdf(3,0,1)))/((covMCS^2)*(1-normcdf(3,0,1)))
            Nsample=(1-(1-normcdf(3,0,1)))/((covMCS^2)*(1-normcdf(3,0,1)));
        end
        for i=1:length(betasegment)-1
            recdummy1_3=sqrt(sum(recdummy1_2.^2,2));
            p1=chi2cdf(betasegment(i)^2,nRV);
            p2=chi2cdf(betasegment(i+1)^2,nRV);
            dummy0=ceil(Nsample*(p2-p1)-length(find(recdummy1_3>betasegment(i) & recdummy1_3<betasegment(i+1))));
            if dummy0<0
                idx=find(recdummy1_3>betasegment(i) & recdummy1_3<betasegment(i+1));
                recdummy1_1(idx(1:abs(dummy0)),:)=[];recdummy1_2(idx(1:abs(dummy0)),:)=[];clear idx
            else
                temp=unifrnd(p1,p2,dummy0,1);
                x=normrnd(0,1,dummy0,nRV);
                x=x.*(sqrt(chi2inv(temp,nRV))./sqrt(sum(x.^2,2)));
                recdummy1_2=vertcat(recdummy1_2,x);recdummy1_1=vertcat(recdummy1_1,NaN(dummy0,1));
            end
        end
        temp=sqrt(sum((recdummy1_2-ucenter).^2,2));
        recdummy1_1(temp<totalbeta(1),:)=[];recdummy1_2(temp<totalbeta(1),:)=[];
        originalusample=recdummy1_2;
        clear x temp p1 p2 dummy0 recdummy1_3 ucenter
        
        %% Subsequent Sphere
        [dummy0,~]=size(candidatepoint1);
        dummy1=0;
        for g=2:dummy0
            if dummy1<nRV
                ucenter=candidatepoint1(1,:);
                candidatepoint1(1,:)=[];
            else
                ucenter=candidatepoint1(1,:);
                candidatepoint1(1,:)=[];
            end
            [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,ucenter);DATATESTING=vertcat(DATATESTING,sample);
            centerFE=G(eco,sample);
            if centerFE<0
                category=0;
            else
                category=1;
            end
            DATASET=vertcat(DATASET,[sample category]);
            FE=FE+1;
            if centerFE<0 || isempty(recfdummy1_2)==1
                continue
            end
            if min(sqrt(sum((recfdummy1_2-ucenter).^2,2)))<1
                continue
            end
            candidatepoint=vertcat(candidatepoint,ucenter);
            
            beta1=min(sqrt(sum((recfdummy1_2-ucenter).^2,2)));beta2=beta1-dring2;betabest=beta1;
            fakerecdummy1_2=recdummy1_2;rectodelete=[];
            while beta1>0
                dummy2=sqrt(sum((fakerecdummy1_2-ucenter).^2,2));
                temp=find(dummy2>beta2 & dummy2<beta1);
                dummy2=dummy2(temp,:);dummy3=fakerecdummy1_2(temp,:);
                if isempty(dummy3)==1
                    break
                end
                dummy6=0;dummy8=[];dummy9=[];[dummy10,~]=size(dummy2);counter=counternum;
                while dummy6<nagent2
                    [~,ia,~]=intersect(temp,rectodelete);
                    roulete=1:1:dummy10;roulete(ia)=[];
                    if isempty(roulete)==1
                        break
                    end
                    dummy4=roulete(ceil(unifrnd(0,length(roulete),1,1)));
                    dummy5=dummy3(dummy4,:);
                    if isnan(recdummy1_1(temp(dummy4)))==0
                        dummy7=recdummy1_1(temp(dummy4));
                        rectodelete=vertcat(rectodelete,temp(dummy4));
                    else
                        if AIstatus==1 && counter<0
                            [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,dummy5);DATATESTING=vertcat(DATATESTING,sample);
                            label=round(net(sample'))';
                            if label==0
                                dummy7=G(eco,sample);
                                if dummy7<0
                                    category=0;
                                else
                                    category=1;
                                end
                                DATASET=vertcat(DATASET,[sample category]);
                                FE=FE+1;
                            else
                                dummy7=1;
                            end
                        else
                            [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,dummy5);DATATESTING=vertcat(DATATESTING,sample);
                            dummy7=G(eco,sample);
                            if dummy7<0
                                category=0;
                            else
                                category=1;
                            end
                            DATASET=vertcat(DATASET,[sample category]);
                            FE=FE+1;
                            counter=counter-1;
                        end
                        recdummy1_1(temp(dummy4))=dummy7;
                        rectodelete=vertcat(rectodelete,temp(dummy4));
                    end
                    if dummy7<0
                        dummy8=vertcat(dummy8,dummy7);dummy9=vertcat(dummy9,dummy5);
                        dummy6=dummy6+1;
                    end
                end
                clear dummy3 dummy4 dummy5 dummy7 temp counter
                
                if dummy6==0
                    break
                end
                
                for h=1:dummy6
                    segmentc(1,:)=0.5;increment=0.5;
                    for i=1:linearfold2
                        increment=increment/2;
                        temp1(i,:)=(dummy9(h,:)-ucenter).*segmentc(i,:)+ucenter;
                        temp2(i,:)=sqrt(sum((temp1(i,:)-ucenter).^2));
                        [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,temp1(i,:));DATATESTING=vertcat(DATATESTING,sample);
                        temp4(i,:)=G(eco,sample);
                        if temp4(i,:)<0
                            category=0;
                        else
                            category=1;
                        end
                        DATASET=vertcat(DATASET,[sample category]);
                        FE=FE+1;
                        if temp4(i,:)<0 && i<linearfold2
                            segmentc(i+1,:)=segmentc(i,:)-increment;
                        else
                            segmentc(i+1,:)=segmentc(i,:)+increment;
                        end
                    end
                    
                    segmentc(end,:)=[];[~,idx]=sort(segmentc);
                    segmentc=segmentc(idx);
                    temp2=temp2(idx);
                    temp4=temp4(idx);
                    
                    %Quadratic FIT
                    temp5=vertcat(centerFE,temp4,dummy8(h));
                    temp6=vertcat(0,temp2,beta1);
                    temp7=vertcat(0,segmentc,1);
                    clear temp1 temp2 temp3 temp4 idx
                    for i=1:length(temp5)
                        if temp5(i)<0
                            break
                        end
                    end
                    p=polyfit(temp6(i-1:i),temp5(i-1:i),1);
                    storedbeta(h)=-p(2)/p(1);
                    p=polyfit(temp7(i-1:i),temp5(i-1:i),1);
                    temp8=-p(2)/p(1);
                    recfdummy1_2=vertcat(recfdummy1_2,(dummy9(h,:)-ucenter).*temp8+ucenter);
                    clear segmentc p temp5 temp6 temp7 temp8
                end
                
                if min(storedbeta)<beta2
                    beta1=min(storedbeta);
                    beta2=beta1-dring2;
                end
            end
            totalbeta=vertcat(totalbeta,beta2);
            clear dummy6 dummy8 dummy9 dummy10 storedbeta beta1 beta2
            dummy1=dummy1+1;
            
            dummy2=sqrt(sum((recdummy1_2-ucenter).^2,2));
            recdummy1_1(dummy2<totalbeta(dummy1+1),:)=[];recdummy1_2(dummy2<totalbeta(dummy1+1),:)=[];
            if dummy1>=sphere
                break
            end
        end
        clear g h i j dummy0 dummy2 fakerecdummy1_2
        
        counter=counternum;
        for i=1:length(recdummy1_1)
            if isnan(recdummy1_1(i))==1
                if AIstatus==1 && counter<0
                    dummy=find(isnan(recdummy1_1)==1);
                    [sample,~,~,~]=summonsample(length(dummy),mu,sigma,nRV,dist,recdummy1_2(dummy,:));DATATESTING=vertcat(DATATESTING,sample);
                    label=round(net(sample'))';
                    recdummy1_1(dummy,:)=label;
                    break
                else
                    [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,recdummy1_2(i,:));DATATESTING=vertcat(DATATESTING,sample);
                    recdummy1_1(i,:)=G(eco,sample);
                    if recdummy1_1(i,:)<0
                        category=0;
                    else
                        category=1;
                    end
                    DATASET=vertcat(DATASET,[sample category]);
                    FE=FE+1;
                    counter=counter-1;
                end
            end
        end
        Pf=length(recdummy1_1(recdummy1_1<=0))/Nsample;
        firstbeta=totalbeta(1);
    end
end

end

