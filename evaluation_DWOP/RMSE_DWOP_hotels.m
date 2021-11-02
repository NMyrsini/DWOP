%Read Data

%tripadvisor_dataset
%Please set the correct file path
fid1 = fopen('tripadvisor_dataset\Preveza_Tripadvisor.txt');
C = textscan(fid1, '%d %d %f %s ', 'delimiter',',');
fclose(fid1);
elements = numel(C{1});
%Load Predictions
data = load('Preveza_Tripadvisor_DWOP_predictions.mat');
vars = whos('-file','Preveza_Tripadvisor_DWOP_predictions.mat');
el = numel(vars);
predicted_ratings = data.PredMean1;

    
real_ratings = (C{3});

RMSE = sqrt(sum((real_ratings-predicted_ratings).^2)/elements);
