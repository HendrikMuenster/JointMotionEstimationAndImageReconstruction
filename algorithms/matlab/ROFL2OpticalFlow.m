function [x,y] = ROFL2OpticalFlow(K,dimsU,tol,sequence,v,alpha,gamma,varargin)
    nPx = prod(dimsU);
    %insert is a list of images, output a list of velocity fields
    %returns primal and dual variable
    
    if (numel(varargin)==1)
        varargin = [varargin{:}];
    end

    v = v(:);
    
    v1 = v(1:nPx);
    v2 = v(nPx+1:end);

    vararginParser;
	
	D3Dp = generateForwardGradientND( dimsU,[1,1,1]);
	D3Dc = generateCentralGradientND( dimsU,[1,1,1]);
	
	Dx3Dc = D3Dc(0*nPx+1:1*nPx,:);
	Dy3Dc = D3Dc(1*nPx+1:2*nPx,:);
	
	Dx3Dp = D3Dp(0*nPx+1:1*nPx,:);
	Dy3Dp = D3Dp(1*nPx+1:2*nPx,:);
	Dt3Dp = D3Dp(2*nPx+1:3*nPx,:);
	
    v1D = spdiags(v1,0,nPx,nPx);
    v2D = spdiags(v2,0,nPx,nPx);
    
    transportPart = Dt3Dp + v1D*Dy3Dc + v2D*Dx3Dc;
	
    spat = [Dx3Dp;Dy3Dp];
	
	clear D3Dp D3Dc Dx3Dc Dy3Dc Dx3Dp Dy3Dp Dt3Dp;
    
    tau=1/(max([1;abs(v)])*sqrt(4*numel(dimsU)));
    sigma = tau;
    
    if (~exist('typeNorm','var'))
        typeNorm = 1;
    end
    
    if (~exist('x','var'))
        x = zeros(prod(dimsU),1);
    end
    if (~exist('y','var'))
        y1 = zeros(size(K,1),1);
        y2 = zeros(size(spat,1),1);
        y3 = zeros(size(transportPart,1),1);
    else
        y1 = y(1:size(K,1));
        y2 = y(size(K,1)+1:size(K,1)+size(spat,1));
        y3 = y(size(K,1)+size(spat,1)+1:end);
    end
    if (~exist('Kx','var'))
        Kx1 = y1;
        Kx2 = y2;
        Kx3 = y3;
    else
        Kx1 = Kx(1:size(K,1));
        Kx2 = Kx(size(K,1)+1:size(K,1)+size(spat,1));
        Kx3 = Kx(size(K,1)+size(spat,1)+1:end);
    end



    %if (~exist('Kty','var'))
        Kty = x;
    %end
    if (~exist('maxInnerIt','var'))
        maxInnerIt = 50000;
    end
    if (~exist('alphaScaling','var'))
        alphaScaling = ones(prod(dimsU),1);
    end

    %rewrite alphaScaling
    alpha = alpha * (alphaScaling(:)+1e-8);
    
    reverseStr = [];
    
    %predefine vector v

    iterations = 1;err=1;
    while err>tol && iterations < maxInnerIt
        y1Old = y1;
        y2Old = y2;
        y3Old = y3;
        
        xOld = x;
        
        if (mod(iterations,100)==1)
            reverseStr = printToCmd( reverseStr,sprintf('ROF-L2 Optical Flow\nIteration: #%d : Residual %f\n',iterations,err) );
        end
        
        KtyOld = Kty;
        
        Kty = (K'*y1 + spat'*y2 + transportPart'*y3);

        %primal prox
        x = x - tau*Kty;
        
        Kx1Old = Kx1;
        Kx2Old = Kx2;
        Kx3Old = Kx3;
        
        Kx1 = K*x;
        Kx2 = spat*x;
        Kx3 = transportPart*x;
        
        %prox operators
        y1 = (y1 + sigma*(2*Kx1-Kx1Old))/(sigma+1) - sigma/(sigma+1)*sequence(:);
		
        if (typeNorm==2)
            y2Tilde = y2 + sigma*(2*Kx2-Kx2Old);
            norm = sqrt(y2Tilde(1:nPx).^2 + y2Tilde(nPx+1:end).^2);
            y2 = y2Tilde ./max(1,[norm./alpha;norm./alpha]);
        else
            y2 = max(-[alpha;alpha],min([alpha;alpha],y2 + sigma*(2*Kx2-Kx2Old)));
        end

        y3 = gamma/(sigma+gamma)* (y3 + sigma*(2*Kx3-Kx3Old)) ;

        xDiff = xOld - x;
        yDiff = [y1Old - y1;y2Old - y2;y3Old - y3];
        
        KxDiff = [Kx1Old-Kx1;Kx2Old-Kx2;Kx3Old-Kx3];

        %primal residul
        p = xDiff/tau - (KtyOld-Kty);
        %dual residual
        d = yDiff/sigma - KxDiff;

        p=sum(abs(p));
        d=sum(abs(d));

        err = (sum(p)+sum(d)) / nPx;

        iterations = iterations + 1;
    end
    
    y = [y1;y2;y3];
	
	x = reshape(x,dimsU);
    
    %clear comand output
    printToCmd( reverseStr,'');
        
end

