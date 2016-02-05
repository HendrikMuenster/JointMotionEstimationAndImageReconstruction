function [x,y] = motionEstimationPyramid(u,dimsU,tol,alpha,algorithmName,numDualVars,varargin)
    nPx = prod(dimsU);

    if (numel(varargin) == 1)
        varargin = [varargin{:}];
    end
    vararginParser;

	initVar('verbose',0);
    initVar('numPrimalVars',2);
    initVar('x',zeros([dimsU,numPrimalVars]));
    initVar('y',zeros([dimsU,numDualVars]));
    initVar('maxIt',10000);
    initVar('typeNorm',3);
    initVar('stepsize',[1 1 1]);
    initVar('discretization',1);
    initVar('numsteps',200);
    initVar('steplength',0.9);
    initVar('numberOfWarps',10);
    
    %option to define useCPP from outside
    if (~exist('useCPP','var'))
        if (exist([algorithmName,'CPP'],'file'))
            useCPP = 1;
        else
            useCPP = 0;
        end
    end

    %automatic step length generation
    steps = 1;
    for i=2:numsteps
        if (steplength*steps(i-1)*size(u,1) > 5 && steplength*steps(i-1)*size(u,2) > 5)
            steps(i) = steplength*steps(i-1);
        else
            break;
        end
    end
    steps = fliplr(steps);
    
    
    strInner = [];

    for i=1:numel(steps)
        if (verbose > 0)
            disp(['Starting zoom-level ',num2str(steps(i))]);
        end

        %no variables set
        if (i == 1)
            oldHeight = size(u,1);
            oldWidth = size(u,2);
        else
            oldHeight = size(uTmp,1);
            oldWidth = size(uTmp,2);
        end

        checkImage = imresize(u(:,:,1),steps(i),'cubic');

        newHeight = size(checkImage,1);
        newWidth = size(checkImage,2);

        newSize = [newHeight newWidth];

        uTmp = [];
        xTmp = [];
        yTmp = [];

        %rescale all variable to current level
        for j=1:dimsU(3)
            uTmp(:,:,j) = imresize(u(:,:,j),newSize,'cubic');

            if (j < dimsU(3))
                %rescale velocity field
                xTmp(:,:,j,1) = imresize(x(:,:,j,1),newSize,'nearest')*newHeight/oldHeight;
                xTmp(:,:,j,2) = imresize(x(:,:,j,2),newSize,'nearest')*newWidth/oldWidth;
                
                for k=3:numPrimalVars
                    xTmp(:,:,j,k) = imresize(x(:,:,j,k),newSize,'nearest');
                end

                for k=1:numDualVars
                    yTmp(:,:,j,k) = imresize(y(:,:,j,k),newSize,'nearest');
                end
            end
        end
        
        x = [];
        y = [];

        for j=1:dimsU(3)-1

            if (useCPP)
                xT = xTmp(:,:,j,:);
                yT = yTmp(:,:,j,:);
            
                if (exist('alpha1','var'))
                    [xRes,yRes] = eval([algorithmName,'CPP','(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,alpha1,maxIt,typeNorm,xT,yT,stepsize,discretization,numberOfWarps)']);
                    
                    xRes = reshape(xRes,[newHeight,newWidth,numPrimalVars]);
                    for k=1:numPrimalVars
                        x(:,:,j,k) = xRes(:,:,k);
                    end
                else
                    if (strcmp(algorithmName,'L2TVBregOpticalFlow'))
                        [x(:,:,j,1),x(:,:,j,2),yRes] = eval(['L2TVBregOpticalFlowCPP','(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,maxIt,numBreg)']);
                    else
                        [x(:,:,j,1),x(:,:,j,2),yRes] = eval([algorithmName,'CPP','(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,maxIt,typeNorm,xT,yT,stepsize,discretization,numberOfWarps)']);
                    end
                    
                    
                end
                
                yRes = reshape(yRes,[newHeight,newWidth,numDualVars]);
                for k=1:numDualVars
                    y(:,:,j,k) = yRes(:,:,k);
                end
            else
                clear xT yT;
                
                for k=1:numPrimalVars
                    xT(:,k) = reshape(xTmp(:,:,j,k),[prod(newSize),1]);
                end
                for k=1:numDualVars
                    yT(:,k) = reshape(yTmp(:,:,j,k),[prod(newSize),1]);
                end
                
                if (exist('alpha1','var'))
                    [xRes,yRes] = eval([algorithmName,'(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,alpha1,''x'',xT,''y'',yT,''maxIt'',maxIt)']);
                else
                    [xRes,yRes] = eval([algorithmName,'(uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,''x'',xT,''y'',yT,''maxIt'',maxIt)']);
                end
                
                for k=1:numPrimalVars
                    x(:,:,j,k) = xRes(:,:,k);
                end
                for k=1:numDualVars
                    y(:,:,j,k) = yRes(:,:,k);
                end

            end
            

        end

        %figure(4);imagesc(colourfulOrientationPlot(x(:,:,1,1),x(:,:,1,2),5));drawnow;

    end
    
    %add zero frames
    x(:,:,end+1,:) = 0;
    y(:,:,end+1,:) = 0;
   
    printToCmd(strInner,'');
end

