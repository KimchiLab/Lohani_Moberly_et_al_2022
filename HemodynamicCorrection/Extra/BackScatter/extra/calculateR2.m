function R2 = calculateR2(ydata,yest)
%calculate coefficient of determination
SSres=sum( (ydata-yest).^2 );
SStot=sum( (ydata-mean(ydata)).^2 );
R2=1-SSres/SStot;

end
