function[r]=ifisherz(z)
%IFISHERZ Inverse Fisher's Z-transform.
% re-transforms Z into the correlation coefficient R.
for i=1:size(z,1)
    for j=1:size(z,2)
        for k=1:size(z,3)
            tmpz=z(i,j,k);
            tmpz=tmpz(:);
            r(i,j,k)=(exp(2*tmpz)-1)./(exp(2*tmpz)+1);
        end
    end
end