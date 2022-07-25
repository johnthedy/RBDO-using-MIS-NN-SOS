function [candidatepoint1,candidatepoint2]=cp(beta,recdummy1_1,recdummy1_2,nRV)
% Candidate Point
candidatepoint=[];
for i=1:nRV
    dummy0=zeros(1,nRV);
    dummy0(i)=beta;
    candidatepoint=vertcat(candidatepoint,dummy0);
end
for i=1:nRV
    dummy0=zeros(1,nRV);
    dummy0(i)=-beta;
    candidatepoint=vertcat(candidatepoint,dummy0);
end
for i=1:1
    dummy0=ones(1,nRV).*beta;
    candidatepoint=vertcat(candidatepoint,dummy0);
end
for i=1:nRV
    dummy1=nchoosek(1:1:nRV,i);
    [a,~]=size(dummy1);
    for j=1:a
        dummy0=ones(1,nRV).*beta;
        dummy0(dummy1(j,:))=-beta;
        candidatepoint=vertcat(candidatepoint,dummy0);
    end
end
[cplength,~]=size(candidatepoint);

% determine candidate point order
cpindex=[];
for i=1:cplength
    dummy0=0;
    rfailureusample=recdummy1_2(recdummy1_1<0,:);
    for j=1:nRV
        dummy0=dummy0+(mean(rfailureusample(:,j))-candidatepoint(i,j)).^2;
    end
    cpindex(i)=(sqrt(dummy0));
end
[~,cporder1]=sort(cpindex,'descend');
[~,cporder2]=sort(cpindex,'ascend');
candidatepoint1=candidatepoint(cporder1,:);
candidatepoint2=candidatepoint(cporder2,:);