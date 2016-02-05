% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.5
% Date: 2015-12-09

% All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and University of Muenster not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
function [u,v] = TVTVOpticalFlowCombined( sequence,alpha,beta,gamma,dimsU,varargin)
    vararginParser;
    
    nPx = prod(dimsU);
	
	initVar('tolU',1e-5);
	initVar('tolV',1e-5);
    initVar('tol',1e-5);
	initVar('maxItOuter',20);
	initVar('isotropeU',3);
	initVar('stepsize',[1 1 1]);
	initVar('maxIterationsU',50000);
	initVar('K',speye(prod(dimsU)));
     
    %create initial u and v
    initialV = zeros([2*nPx,1]);
    if (exist('ROFL2OpticalFlowCPP','file'))
        [u,yU] = ROFL2OpticalFlowCPP(sequence,tolU,alpha,gamma,maxIterationsU,initialV,stepsize,isotropeU,char);
    else
        [u,yU] = ROFL2OpticalFlow(K,dimsU,tolU,sequence,initialV,alpha,gamma);
    end
    [v,yV] = motionEstimationPyramid(u,dimsU,tolV,beta/gamma,'L1TVOpticalFlow',4,'maxIt',1000);
    
    reverseStr = [];
    
    iterations=1;err = Inf;
    while err > tol && iterations < maxItOuter
		%figure(1001);imagesc(u(:,:,1))
		%figure(1002);imagesc(u(:,:,2))
        %figure(1003);imagesc(u(:,:,3))
        
		%figure(1004);imagesc(flowToColorV2(cat(3,v(:,:,1,1),v(:,:,1,2))))
        %figure(1005);imagesc(flowToColorV2(cat(3,v(:,:,2,1),v(:,:,2,2))))
        %drawnow
        
        uOld = u;
        vOld = v;
        
        reverseStr = printToCmd( [],sprintf('TV-TV-L1 Optical Flow\nIteration: #%d : Residual %f\n',iterations,err) );
        
        %use zooming
        TmpVarargin = varargin;
        TmpVarargin{end+1} = 'x';
        TmpVarargin{end+1} = v;
        TmpVarargin{end+1} = 'y';
        TmpVarargin{end+1} = yV;
        TmpVarargin{end+1} = 'stepsize';
        TmpVarargin{end+1} = stepsize;
        TmpVarargin{end+1} = 'discretization';
        TmpVarargin{end+1} = 1;
        TmpVarargin{end+1} = 'maxIt';
        TmpVarargin{end+1} = 500;
        TmpVarargin{end+1} = 'numsteps';
        TmpVarargin{end+1} = 1;
        
        [v,yV] = motionEstimationPyramid(u,dimsU,tolV,beta/gamma,'L1TVOpticalFlow',4,TmpVarargin);
        
        %check for CPP version
        if (exist('ROFL1OpticalFlowCPP','file'))
            [u,yU] = ROFL1OpticalFlowCPP(sequence,tolU,alpha,gamma,maxIterationsU,v,stepsize,isotropeU,ones(dimsU),u,yU);
        else
            [u,yU] = ROFL1OpticalFlow(K,dimsU,tolU,sequence,v,alpha,gamma,'x',reshape(u,[nPx,1]),'y',yU);
        end
        
        err = (sum(abs(u(:)-uOld(:))) / nPx + sum(abs(v(:)-vOld(:))) / (2*nPx)) / 2;
        
        iterations = iterations + 1;
	end
    
    printToCmd( reverseStr,'');

     u = reshape(u,dimsU);
     v = reshape(v,[dimsU,2]);
    
end

