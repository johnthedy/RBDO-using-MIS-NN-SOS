function [sample,usample,dummy1,dummy2]=summonsample(n,mu,sigma,nRV,dist,u)
for i=1:length(mu)
    if strcmp(dist{i},'Normal')==1
        dummy1(i)=mu(i);
        dummy2(i)=sigma(i);
    elseif strcmp(dist{i},'Lognormal')==1
        dummy2(i)=sqrt(log(1+(sigma(i)^2)/(mu(i)^2)));
        dummy1(i)=log(mu(i))-0.5*dummy2(i)^2;
    elseif strcmp(dist{i},'Extreme Value')==1
        dummy2(i)=sigma(i)/1.2825;
        dummy1(i)=(mu(i)+dummy2(i)*0.5772);
    elseif strcmp(dist{i},'Gamma')==1
        dummy1(i)=(mu(i)^2)/(sigma(i)^2);
        dummy2(i)=(sigma(i)^2)/mu(i);
    end
end

%Sample of random variable
sample=zeros(n,nRV);usample=zeros(n,nRV);
for i = 1:n
    for j=1:nRV
        dummy0=cdf('Normal',u(i,j),0,1); %Random number from 0-1
        if dummy0==1
            dummy0=1-1e-10;
        end
        sample(i,j)=icdf(dist{j},dummy0,dummy1(j),dummy2(j)); %Generate random sample in original space based on earlier random number
        usample(i,j)=icdf('Normal',dummy0,0,1); %Generate random sample in standard normal space based on earliner random number
    end
end
end

