function Population = EnvironmentalSelection(Population,V,theta,raw)

PopObj = Population.objs;
[N,M]  = size(PopObj);
NV     = size(V,1);

%% Translate the population
PopObj = PopObj - repmat(min(PopObj,[],1),N,1);

%% Calculate the degree of violation of each solution
if raw
    CV = sum(max(0,Population.cons),2);
else
    CV = zeros(N,1);
end

%% Calculate the smallest angle value between each vector and others
cosine = 1 - pdist2(V,V,'cosine');
cosine(logical(eye(length(cosine)))) = 0;
gamma  = min(acos(cosine),[],2);

%% Associate each solution to a reference vector
Angle = acos(1-pdist2(PopObj,V,'cosine'));
[~,associate] = min(Angle,[],2);

%% Select one solution for each reference vector
Next = zeros(1,NV);
next = [];
for i = unique(associate)'
    current = find(associate==i);
    [FrontNo,MaxFNo] = NDSort(PopObj(current,:),CV(current,:),N);
    tmp = current(find(FrontNo == 1));
    if ~isempty(tmp)
        % Calculate the APD value of each solution
        APD = (1+M*theta*Angle(tmp,i)/gamma(i)).*sqrt(sum(PopObj(tmp,:).^2,2));
        % Select the one with the minimum APD value
        [~,best] = min(APD);
    end
    next = [next,tmp(best)];
    
end
% Population for next generation
Population = Population(next);
end