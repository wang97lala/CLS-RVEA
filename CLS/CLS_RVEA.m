classdef CLS_RVEA < ALGORITHM
% <multi/many> <real> <large/none> <constrained/none>
% L
    methods
        function main(Algorithm,Problem)
            
            
            %% Generate random population
            k = Algorithm.ParameterSet(5);
            [V,Problem.N] = UniformPoint(Problem.N,Problem.M);

            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();

            Population1    = EnvironmentalSelection(Population1,V,(Problem.FE/Problem.maxFE)^2,true);
            Population2    = EnvironmentalSelection(Population2,V,(Problem.FE/Problem.maxFE)^2,false);
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                Fitness = calFitness(Population1.objs);
                if length(Population1) >= 2
                    Rank = randperm(length(Population1),floor(length(Population1)/2)*2);
                else
                    Rank = [1,1];
                end
                Loser  = Rank(1:end/2);
                Winner = Rank(end/2+1:end);
                Change = Fitness(Loser) >= Fitness(Winner);
                Temp   = Winner(Change);
                Winner(Change) = Loser(Change);
                Loser(Change)  = Temp;
                Offspring1      = Operator(Population1(Loser),Population1(Winner));
                
                
                Fitness = calFitness(Population2.objs);
                if length(Population2) >= 2
                    Rank = randperm(length(Population2),floor(length(Population2)/2)*2);
                else
                    Rank = [1,1];
                end
                Loser  = Rank(1:end/2);
                Winner = Rank(end/2+1:end);
                Change = Fitness(Loser) >= Fitness(Winner);
                Temp   = Winner(Change);
                Winner(Change) = Loser(Change);
                Loser(Change)  = Temp;
                Offspring2      = Operator(Population2(Loser),Population2(Winner));
                
                
                Population1     = EnvironmentalSelection([Population1,Offspring1,Offspring2],V,(Problem.FE/Problem.maxFE)^2,true);
                Population2     = EnvironmentalSelection([Population2,Offspring2,Offspring1],V,(Problem.FE/Problem.maxFE)^2,false);
                
            end
        end
    end
end

function Fitness = calFitness(PopObj)
% Calculate the fitness by shift-based density

    N      = size(PopObj,1);
    fmax   = max(PopObj,[],1);
    fmin   = min(PopObj,[],1);
    PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    Dis    = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    Fitness = min(Dis,[],2);
end