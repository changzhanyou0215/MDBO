

function [fMin, bestX, Convergence_curve] = MDBO(pop, M, c, d, dim, fobj)

P_percent = 0.2; % Proportion of dung beetles performing ball-rolling behavior  
pNum = round(pop * P_percent); % Number of dung beetles performing ball-rolling behavior 
lb = c.*ones(1,dim);    % lower bound vector for the value
ub = d.*ones(1,dim);    % upper bound vector for the value

% ★★ Improvement 1: Initialization of ROBL populations ★★
Sol_forward=initialization_for(pop,dim,ub,lb);
Sol_backward=initialization_back(Sol_forward,pop,dim,ub,lb);
Sol_all=[Sol_forward;Sol_backward];
for i = 1:2*pop
    Sol_all_fitness(i)=fobj(Sol_all(i,:));
end
[~,sorted_indexes]=sort(Sol_all_fitness);
for i = 1:pop
    x(i,:)=Sol_all(sorted_indexes(i),:);
    fit(i)=fobj(x(i,:));
end


pFit = fit;                       
pX = x; 
XX = pX;    
[fMin, bestI] = min(fit); % fMin is the value of global optimal fitness
bestX = x(bestI,:);       % bestX is the global optimal solution corresponding to fMin

% Beginning iteration for individual updates
for t = 1 : M 

[fmax, B] = max(fit);
worse = x(B,:);   
r2 = rand(1);

% Rolling Ball Behavior Dung Beetle Position Updates (divided into No Obstacle Mode and Obstacle Mode)
for i = 1 : pNum    
    if(r2<0.9)
        % ① No Obstacle Mode
        % ★★ Improvement 2: Gold Sine Strategy ★★
        R1 = rand()*2*pi;
        R2 = rand()*pi;
        tao = (sqrt(5)-1)/2;
        x1 = -pi+(1-tao)*2*pi;
        x2 = -pi+tao*2*pi; 
        x(i,:) = pX(i,:)*abs(sin(R1))+R2*sin(R1)*abs(x1*pX(i,:)-x2*x(i,:));
    else
        % ② Obstacle Mode
        aaa = randperm(180,1);
        if (aaa==0 ||aaa==90 ||aaa==180)
            x(i,:) = pX(i,:);   
        end
        theta = aaa*pi/180;   
        x(i,:) = pX(i,:)+tan(theta).*abs(pX(i,:)-XX(i,:)); 
    end
    % 越界校正
    x(i,:) = Bounds(x(i,:), lb, ub);    
    fit(i) = fobj(x(i,:));
end

[fMMin, bestII] = min(fit);   % fMin is the current optimal fitness value
bestXX = x(bestII,:);         % bestXX is the current optimum solution
R = 1-t/M;                         

Xnew1 = bestXX.*(1-R); 
Xnew2 = bestXX.*(1+R);                    
Xnew1 = Bounds(Xnew1, lb, ub);
Xnew2 = Bounds(Xnew2, lb, ub);

% Reproductive Behavior of Dung Beetles Location Updates
for i = (pNum + 1) : 12
    x(i,:) = bestXX+((rand(1,dim)).*(pX(i,:)-Xnew1)+(rand(1,dim)).*(pX(i,:)-Xnew2));
    x(i,:) = Bounds(x(i,:), Xnew1, Xnew2);
    fit(i) = fobj(x(i,:)) ;
end

Xnew11 = bestX.*(1-R); 
Xnew22 = bestX.*(1+R);                    
Xnew11 = Bounds(Xnew11, lb, ub);
Xnew22 = Bounds(Xnew22, lb, ub);

% Foraging Behavior of Small Dung Beetles Location Updates
for i = 13 : 19        
    x(i,:) = pX(i,:)+((randn(1)).*(pX(i,:)-Xnew11)+((rand(1,dim)).*(pX(i,:)-Xnew22)));
    x(i,:) = Bounds(x(i,:), lb, ub);
    fit(i) = fobj(x(i,:));
end

% Stealing Behavior Dung Beetle Location Update
for i = 20 : pop
    x(i,:) = bestX + randn(1,dim).*((abs((pX(i,:)-bestXX)))+(abs((pX(i,:)-bestX))))./2;
%     % 引入模拟退火算法
% 	x(i,:) =sim_anl(fobj, x(i,:), lb, ub, 10, 10E-4);
    x(i,:) = Bounds(x(i,:), lb, ub);
    fit(i) = fobj(x(i,:)) ;
end

%% 4. Improvement point: introduction of t-distribution variation

for i = 1:pop
    A = (std(fit) - mean(fit))/(mean(fit)^2);
    if(t<=M/2 && A>10E-2)
        % ori_value = rand(1,dim);
        % cauchy_value = tan((ori_value-0.5)*pi);
        Temp = x( i, : ) + x( i, : )*trnd(i); % Variation of t-distribution based on number of iterations
        % Temp = x( i, : ) + x( i, : ).*cauchy_value; % Introduction of Cauchy variation
        % border processing
        Temp(Temp>ub) = ub(Temp>ub);
        Temp(Temp<lb) = lb(Temp<lb);
        fitvalue = fobj(Temp);
        if(fitvalue <fit(i))
            x( i, : ) = Temp;
            fit(i) = fitvalue;
        end
    end
end

% Introducing a simulated annealing algorithm
for i = 1 : pop
    % Introducing a simulated annealing algorithm
	x(i,:) =sim_anl(fobj, x(i,:), lb, ub, 10, 10E-4);
    x(i,:) = Bounds(x(i,:), lb, ub);
    fit(i) = fobj(x(i,:)) ;
end


% Update the individual optimum and the global optimum
XX = pX;
for i = 1 : pop 
    if (fit(i) < pFit(i))
        pFit(i) = fit(i);
        pX(i,:) = x(i,:);
    end

    if( pFit(i) < fMin)
        fMin = pFit(i);
        bestX = pX(i,:);
    end
end

% Iteration Curve Records
Convergence_curve(t) = fMin;
end

% Out-of-bounds correction function
function s = Bounds(s, Lb, Ub)
temp = s;
I = temp < Lb;
temp(I) = Lb(I);

J = temp > Ub;
temp(J) = Ub(J);
s = temp;

