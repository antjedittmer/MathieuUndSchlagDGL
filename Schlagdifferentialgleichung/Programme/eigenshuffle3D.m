function [eigval,eigvec] = eigenshuffle3D(A)


r=size(A);
if (r(1)~=r(2))
   error('The argument must be a 3D array pxpxn');
   %else
end

if (length(r)~=3)
   %error('The argument must be a 3D array pxpxn');
   %else
   r(3) = 1;
end


if (length(r)~=3) || (r(1)~=r(2))
    error('The argument must be a 3D array pxpxn');
else
    eigvec=zeros(r(1),r(2),r(3));
    eigval=zeros(r(3),r(1));
    for k=1:r(3)
        [evc,evl]=eig(A(:,:,k));
        evl=diag(evl)';
        if k~=1
            dotp=abs(eigvec(:,:,k-1)'*evc); % Multiplizier EVec der vorherigen mit der jetzigen A-Matrix: Evec2 = Evec(k-1)Evec(k)
            [r,idx]=max(dotp); % Maximum Wert jedes Evec2 pro Spalte: idx: Index vom Maximum
            for i=1:length(idx)
                doppi=find(i==idx,2); % Doppelte 'Maximal'Werte
                if length(doppi)==2 %'Handling von
                    for j=1:length(idx)
                        r=find(j==idx,1);
                        if isempty(r)
                            idx(doppi(1))=j;
                        end
                    end
                    for j=k-1:-1:1
                        if eigval(j,idx(doppi(1)))~=eigval(j,idx(doppi(2)))
                            break;
                        end
                    end
                    [r,ind1]=max(evl(doppi));
                    [r,ind2]=max(eigval(j,idx(doppi)));
                    if ind1~=ind2
                        idx(doppi)=fliplr(idx(doppi));
                    end
                end
            end
            for i=1:length(evl)
                eigval(k,idx(i))=evl(i);
                eigvec(:,idx(i),k)=evc(:,i);
            end
        else
            eigvec(:,:,k)=evc;
            eigval(k,:)=evl;
        end
    end
end