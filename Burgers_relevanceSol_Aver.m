clear 
clc 
Tem=[10,20,40,80,160,320];

%[-1,1]上的高斯点和权重

RelevantGaussianPoints = [-0.978228658146057, -0.887062599768095, -0.730152005574049, -0.519096129206812, -0.269543155952345, 0, ...
                   0.269543155952345, 0.519096129206812, 0.730152005574049, 0.887062599768095, 0.978228658146057];

weights = [0.055668567116174, 0.125580369464905, 0.186290210927734, 0.233193764591990, 0.262804544510247, 0.272925086777901, ...
           0.262804544510247, 0.233193764591990, 0.186290210927734, 0.125580369464905, 0.055668567116174];    


% 定义区间数和高斯点数
numIntervals      = Tem(1);
numGaussianPoints = length(RelevantGaussianPoints );


% 生成[0, 2]范围内的等间隔向量
x = linspace(0, 2, numIntervals + 1);

% 计算每个区间的高斯点
gaussianPoints     = zeros(numIntervals, numGaussianPoints);
gaussianValues     = zeros(numIntervals, numGaussianPoints);
gaussianQuadrature = zeros(numIntervals, 1);
cellAverage        = zeros(numIntervals, 1);

for i = 1:numIntervals
    gaussianPoints(i, :) = (x(i + 1) + x(i)) / 2  + (x(i + 1) - x(i)) / 2 * RelevantGaussianPoints;
end

for i = 1 : numIntervals
    for j = 1 : numGaussianPoints
        %非线性方程的隐函数解
        equation = @(y) y - sin(pi *gaussianPoints(i,j)- y / 2)-1/2;

        % 使用fsolve求解方程
        options = optimoptions(@fsolve,'OptimalityTolerance',1.000e-10,'FunctionTolerance',1.000e-10);
        gaussianValues(i,j) = fsolve(equation,0,options );
    end
end

for i = 1 : numIntervals
    sum = 0;
    for j = 1 : numGaussianPoints
        sum = sum + weights(j) * gaussianValues(i,j);
    end
    gaussianQuadrature(i,:) = (x(i + 1) - x(i)) / 2 * sum;
    cellAverage(i,:)        = gaussianQuadrature(i,:) / (x(i + 1) - x(i));
end
cellAverage = cellAverage';

disp('区间的单元平均：');
disp(cellAverage);


 