function CVs = gridConvergenceTester(input1, input2, operator_split)

grid_counts = [10, 20, 30, 40, 50, 75, 100, 125, 150, 200, 500, 1000];

CVs = zeros(length(grid_counts),1);
for k = 1:length(grid_counts)
    
    createBlankTestProblem('test',grid_counts(k));
    pause(2);
    load('test.mat','problem');
    if operator_split
        CVs(k) = calculateCV_OperatorSplit(problem, 'planar', input1, input2);
    else
        CVs(k) = calculateCV(problem, 'planar', input1, input2);
    end
    fprintf('Completed run with %d gridpoints\n', grid_counts(k));
    
end

