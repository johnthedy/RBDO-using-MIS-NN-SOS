%SOS1 means data is accumulated, the result is named Result 1
clc;clear;tic;
restoredefaultpath
%% Opt parameter   
ite=50;
ecosize=10;

%% Random and Optimized Variable
muorigin=[2500 250 125 40];
sigmaorigin=muorigin.*[0.2 0.3 0.3 0.1];
covparam=0.1;
dist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal'};
nRV=6;nopt=2;
betatarget=3;
ub=ones(1,nopt)*2;
lb=ones(1,nopt)*0.01;
Rub=sqrt(chi2inv(1-(1-normcdf(betatarget,0,1)),nRV));
Rlb=0;

%% Reliability Parameter
dring1=0.1;
dring2=0.2;
nagent1_1=nRV*5;
nagent1_2=1;
nagent2=1;
stopval=4;
tol1=3.1;
tol2=2.9;
covMCS1=0.2;
covMCS2=0.05;
lowlimitfactor=10;

%% InitialNN
LayerSize=10;
trainsize=10;initialFE=0;
net=patternnet(LayerSize);
net.trainParam.showWindow=0;
net.divideparam.trainRatio=1;
net.divideparam.valRatio=0;
net.divideparam.testRatio=0;

for i=1:trainsize
    eco=rand(1,nopt).*(ub-lb)+lb;
    x=normrnd(0,1,1,nRV);
    mu=horzcat(muorigin,eco);sigma=horzcat(sigmaorigin,eco.*covparam);
    [sample,usample,~,~]=summonsample(1,mu,sigma,nRV,dist,x);
    dummy1=G(eco,sample);initialFE=initialFE+1;
    if dummy1<0
        category=0;
    else
        category=1;
    end
    DATASET(i,:)=[sample category];
end
[net,~]=train(net,DATASET(:,1:nRV)',DATASET(:,nRV+1)');

%% Start first iteration of SOS
eco=zeros(ecosize,nopt);
fitness=zeros(ecosize,1);beta=zeros(ecosize,1);phi=zeros(ecosize,1);rectransition=0;FE=0;
AIstatus=1;
covMCS=covMCS1;

for i=1:ecosize
    eco(i,:)=rand(1,nopt).*(ub-lb)+lb;
    fitness(i,:)=fobj1(eco(i,:));
    %Setup svm model
    mu=horzcat(muorigin,eco(i,:));sigma=horzcat(sigmaorigin,eco(i,:).*covparam);
    [Pf,dummy1,firstbeta,dummy2,~]=fobj2(eco(i,:),dring1,dring2,nagent1_1,nagent1_2,covMCS,nagent2,stopval,mu,sigma,nRV,dist,Rub,Rlb,AIstatus,net);
    beta(i,:)=-norminv(Pf,0,1);
    originbeta(i,:)=firstbeta;
    if abs(Pf-(1-normcdf(betatarget,0,1)))<lowlimitfactor*(1-normcdf(betatarget,0,1)) || Pf<(1-normcdf(betatarget,0,1))
        phi(i,:)=abs(Pf-(1-normcdf(betatarget,0,1)));
    else
        phi(i,:)=-inf;
    end
    FE=FE+dummy1;
    DATASET=vertcat(DATASET,dummy2);
    
    clear dummy1 dummy2
end
[net,~]=train(net,DATASET(:,1:nRV)',DATASET(:,nRV+1)');

%%
counter=0;
for h=1:ite
    DATASET=[];
    DATATESTING=[];a
    counter=counter+1;
    if counter>3 || h==ite
        disp('start calibrating')
        for j=1:ecosize
            fitness(j,:)=fobj1(eco(j,:));
            %Setup svm model
            mu=horzcat(muorigin,eco(j,:));sigma=horzcat(sigmaorigin,eco(j,:).*covparam);
            [Pf,dummy1,firstbeta,dummy2,dummy3]=fobj2(eco(j,:),dring1,dring2,nagent1_1,nagent1_2,covMCS,nagent2,stopval,mu,sigma,nRV,dist,Rub,Rlb,AIstatus,net);
            beta(j,:)=-norminv(Pf,0,1);
            originbeta(j,:)=firstbeta;
            if abs(Pf-(1-normcdf(betatarget,0,1)))<lowlimitfactor*(1-normcdf(betatarget,0,1)) || Pf<(1-normcdf(betatarget,0,1))
                phi(j,:)=abs(Pf-(1-normcdf(betatarget,0,1)));
            else
                phi(j,:)=-inf;
            end
            FE=FE+dummy1;
            DATATESTING=vertcat(DATATESTING,dummy3);
            DATASET=vertcat(DATASET,dummy2);
            
            clear dummy1 dummy2 dummy3
        end
        counter=0;
        disp('finish calibrating')
    end
    
    for i=1:ecosize
        % Update the best Organism
        % find gbest particle
        if sum(beta>=betatarget)>0
            feasible_particle = find(beta>=betatarget);
            [~,best_particle]=min(fitness(feasible_particle));
            best_particle=feasible_particle(best_particle);
        else
            infeasible_particle=find(abs(phi)==min(abs(phi)));
            [~,best_particle]=min(fitness(infeasible_particle));
            best_particle=infeasible_particle(best_particle);
        end
        
        %check if it fulfilled requirement to shift into stage 2
        if length(find(beta<tol1 & beta>tol2))>ecosize/1.5
            if covMCS==covMCS1
                rectransition=h;
            end
            covMCS=covMCS2;
        end
        bestFitness=fitness(best_particle);bestOrganism=eco(best_particle,:);bestbeta=beta(best_particle);
        
        %% Mutualism Phase
        j=i;
        while i==j
            seed=randperm(ecosize); 
            j=seed(1);                  
        end
        % Determine Mutual Vector & Beneficial Factor
        mutualVector=mean([eco(i,:);eco(j,:)]);
        BF1=round(1+rand); BF2=round(1+rand);
        % Calculate new solution after Mutualism Phase
        ecoNew1=eco(i,:)+rand(1,nopt).*(bestOrganism-BF1.*mutualVector); 
        ecoNew2=eco(j,:)+rand(1,nopt).*(bestOrganism-BF2.*mutualVector);
        ecoNew1=bound(ecoNew1,ub,lb); 
        ecoNew2=bound(ecoNew2,ub,lb);
        % Evaluate the fitness of the new solution
        fitnessNew1=fobj1(ecoNew1);
        %Setup svm model
        mu=horzcat(muorigin,ecoNew1);sigma=horzcat(sigmaorigin,ecoNew1.*covparam);
        [Pf,dummy1,firstbeta1,dummy2,dummy3]=fobj2(ecoNew1,dring1,dring2,nagent1_1,nagent1_2,covMCS,nagent2,stopval,mu,sigma,nRV,dist,Rub,Rlb,AIstatus,net);
        betaNew1=-norminv(Pf,0,1);
        if abs(Pf-(1-normcdf(betatarget,0,1)))<lowlimitfactor*(1-normcdf(betatarget,0,1)) || Pf<(1-normcdf(betatarget,0,1))
            phiNew1=abs(Pf-(1-normcdf(betatarget,0,1)));
        else
            phiNew1=-inf;
        end
        FE=FE+dummy1;
        DATATESTING=vertcat(DATATESTING,dummy3);
        DATASET=vertcat(DATASET,dummy2);
        
        clear dummy1 dummy2 dummy3
        fitnessNew2=fobj1(ecoNew2);
        %Setup svm model
        mu=horzcat(muorigin,ecoNew2);sigma=horzcat(sigmaorigin,ecoNew2.*covparam);
        [Pf,dummy1,firstbeta2,dummy2,dummy3]=fobj2(ecoNew2,dring1,dring2,nagent1_1,nagent1_2,covMCS,nagent2,stopval,mu,sigma,nRV,dist,Rub,Rlb,AIstatus,net);
        betaNew2=-norminv(Pf,0,1);
        if abs(Pf-(1-normcdf(betatarget,0,1)))<lowlimitfactor*(1-normcdf(betatarget,0,1)) ||  Pf<(1-normcdf(betatarget,0,1))
            phiNew2=abs(Pf-(1-normcdf(betatarget,0,1)));
        else
            phiNew2=-inf;
        end
        FE=FE+dummy1;
        DATATESTING=vertcat(DATATESTING,dummy3);
        DATASET=vertcat(DATASET,dummy2);
        
        clear dummy1 dummy2 dummy3
        
        % Accept the new solution if the fitness is better
        if fitnessNew1<fitness(i) && betaNew1>=betatarget
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
            phi(i)=phiNew1;
            beta(i)=betaNew1;
            originbeta(i)=firstbeta1;
        elseif abs(phiNew1)<abs(phi(i)) && covMCS==covMCS1
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
            phi(i)=phiNew1;
            beta(i)=betaNew1;
            originbeta(i)=firstbeta1;
        elseif abs(phiNew1)<abs(phi(i)) && covMCS==covMCS2 && beta(i)<betatarget
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
            phi(i)=phiNew1;
            beta(i)=betaNew1;
            originbeta(i)=firstbeta1;
        end
        
        if fitnessNew2<fitness(j) && betaNew2>=betatarget
            fitness(j)=fitnessNew2;
            eco(j,:)=ecoNew2;
            phi(j)=phiNew2;
            beta(j)=betaNew2;
            originbeta(j)=firstbeta2;
        elseif abs(phiNew2)<abs(phi(j)) && covMCS==covMCS1
            fitness(j)=fitnessNew2;
            eco(j,:)=ecoNew2;
            phi(j)=phiNew2;
            beta(j)=betaNew2;
            originbeta(j)=firstbeta2;
        elseif abs(phiNew2)<abs(phi(j)) && covMCS==covMCS2 && beta(j)<betatarget
            fitness(j)=fitnessNew2;
            eco(j,:)=ecoNew2;
            phi(j)=phiNew2;
            beta(j)=betaNew2;
            originbeta(j)=firstbeta2;
        end
        % End of Mutualism Phase 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %% Commensialism Phase
        j=i;
        while i==j
            seed=randperm(ecosize); 
            j=seed(1);                  
        end
        % Calculate new solution after Commensalism Phase    
        ecoNew1=eco(i,:)+(rand(1,nopt)*2-1).*(bestOrganism-eco(j,:));
        ecoNew1=bound(ecoNew1,ub,lb);
        % Evaluate the fitness of the new solution
        fitnessNew1=fobj1(ecoNew1);
        mu=horzcat(muorigin,ecoNew1);sigma=horzcat(sigmaorigin,ecoNew1.*covparam);
        [Pf,dummy1,firstbeta1,dummy2,dummy3]=fobj2(ecoNew1,dring1,dring2,nagent1_1,nagent1_2,covMCS,nagent2,stopval,mu,sigma,nRV,dist,Rub,Rlb,AIstatus,net);
        betaNew1=-norminv(Pf,0,1);
        if abs(Pf-(1-normcdf(betatarget,0,1)))<lowlimitfactor*(1-normcdf(betatarget,0,1)) || Pf<(1-normcdf(betatarget,0,1))
            phiNew1=abs(Pf-(1-normcdf(betatarget,0,1)));
        else
            phiNew1=-inf;
        end
        FE=FE+dummy1;
        DATATESTING=vertcat(DATATESTING,dummy3);
        DATASET=vertcat(DATASET,dummy2);
        
        clear dummy1 dummy2 dummy3
        % Accept the new solution if the fitness is better
        if fitnessNew1<fitness(i) && betaNew1>=betatarget
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
            phi(i)=phiNew1;
            beta(i)=betaNew1;
            originbeta(i)=firstbeta1;
        elseif abs(phiNew1)<abs(phi(i)) && covMCS==covMCS1
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
            phi(i)=phiNew1;
            beta(i)=betaNew1;
            originbeta(i)=firstbeta1;
        elseif abs(phiNew1)<abs(phi(i)) && covMCS==covMCS2 && beta(i)<betatarget
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
            phi(i)=phiNew1;
            beta(i)=betaNew1;
            originbeta(i)=firstbeta1;
        end
        % End of Commensalism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Parasitism Phase
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        % Determine Parasite Vector & Calculate the fitness
        parasiteVector=eco(i,:);
        seed=randperm(nopt);           
        pick=seed(1:ceil(rand*nopt));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(ub(pick)-lb(pick))+lb(pick);
        fitnessParasite=fobj1(parasiteVector);
        %Setup svm model
        mu=horzcat(muorigin,parasiteVector);sigma=horzcat(sigmaorigin,parasiteVector.*covparam);
        [Pf,dummy1,firstbetaparasite,dummy2,dummy3]=fobj2(parasiteVector,dring1,dring2,nagent1_1,nagent1_2,covMCS,nagent2,stopval,mu,sigma,nRV,dist,Rub,Rlb,AIstatus,net);
        betaNew1=-norminv(Pf,0,1);
        if abs(Pf-(1-normcdf(betatarget,0,1)))<lowlimitfactor*(1-normcdf(betatarget,0,1))
            phiNew1=abs(Pf-(1-normcdf(betatarget,0,1)));
        elseif Pf<(1-normcdf(betatarget,0,1))
            phiNew1=(1-normcdf(betatarget,0,1))-Pf;
        else
            phiNew1=-inf;
        end
        FE=FE+dummy1;
        DATATESTING=vertcat(DATATESTING,dummy3);
        DATASET=vertcat(DATASET,dummy2);
        
        clear dummy1 dummy2 dummy3

        % Kill organism j and replace it with the parasite 
        % if the fitness is lower than the parasite
        if fitnessParasite<fitness(j) && betaNew1>=betatarget
            fitness(j)=fitnessParasite;
            eco(j,:)=parasiteVector;
            phi(j)=phiNew1;
            beta(j)=betaNew1;
            originbeta(j)=firstbetaparasite;
        elseif abs(phiNew1)<abs(phi(j)) && covMCS==covMCS1
            fitness(j)=fitnessParasite;
            eco(j,:)=parasiteVector;
            phi(j)=phiNew1;
            beta(j)=betaNew1;
            originbeta(j)=firstbetaparasite;
        elseif abs(phiNew1)<abs(phi(j)) && covMCS==covMCS2 && beta(j)<betatarget
            fitness(j)=fitnessParasite;
            eco(j,:)=parasiteVector;
            phi(j)=phiNew1;
            beta(j)=betaNew1;
            originbeta(j)=firstbetaparasite;
        end
        % End of Parasitism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    tic;[net,~]=train(net,DATASET(:,1:nRV)',DATASET(:,nRV+1)');
    recNNTtime(h,:)=toc;
    recNNdata(h,:)=length(DATASET(:,1));
    recDATATESTING{h}=DATATESTING(:,5:6);
    fprintf('COV of MCS = %d\nTotal Performed Function Evaluation Up To Now = %d\n',covMCS,FE)
    disp('beta, firstsphere radius, difference toward Pf target, and fitness :')
    disp([beta originbeta phi fitness])
    
    recbestfitness(h)=bestFitness;
    recbestbeta(h)=bestbeta;
    recFE(h)=FE;
    receco{h}=eco;
    
    if h==1
        recFEperRA(h)=FE/(4*ecosize+ecosize);
    elseif h==rectransition
        recFEperRA(h)=(FE-recFE(h-1))/(4*ecosize+ecosize);
    else
        recFEperRA(h)=(FE-recFE(h-1))/(4*ecosize);
    end
    
end
fprintf('bestFitness = %d',bestFitness)
fprintf('bestOrganism = %d',bestOrganism)
