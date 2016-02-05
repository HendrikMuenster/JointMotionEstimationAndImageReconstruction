clear all;close all;clc;

addpath(genpath(cd));

load('data/data.mat');

gtFlow = readFlowFile('data/stuff_rof.flo');

for i=1:size(imageSequence,3)
    figure(1);imagesc(imageSequence(:,:,i));title('Ground-truth image')
    figure(2);imagesc(imageSequenceNoise(:,:,i));title('Noisy image')
    figure(3);imagesc(flowToColorV2(cat(3,gtFlow(:,:,1),gtFlow(:,:,2))));title('Ground-truth flow')
end

%%

alpha = 0.015; %regularization weight for images
beta = 0.05; %regularization weight for velocity field
gamma = 1; %weight for optical flow constraint

[u,v] = TVTVOpticalFlowCombined(imageSequenceNoise,alpha,beta*gamma,gamma,size(imageSequence),'maxItOuter',100);

%%
for i=1:size(imageSequence,3)
    figure(1);imagesc(imageSequence(:,:,i));title('Ground-truth image')
    figure(2);imagesc(imageSequenceNoise(:,:,i));title('Noisy image')
    figure(3);imagesc(u(:,:,i));title('Reconstructed image')
    figure(4);imagesc(flowToColorV2(cat(3,gtFlow(:,:,1),gtFlow(:,:,2))));title('Ground-truth flow')
    figure(5);imagesc(flowToColorV2(cat(3,v(:,:,i,1),v(:,:,i,2))));title('Reconstructed flow')
    
    pause
end