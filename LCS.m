function [D, dist, aLongestString] = LCS(X,Y)

n =length(X);
m =length(Y);
L=zeros(n+1,m+1);
L(1,:)=0;
L(:,1)=0;
b = zeros(n+1,m+1);
b(:,1)=1;%%%Up
b(1,:)=2;%%%Left

for i = 2:n+1
    for j = 2:m+1
        if (X(i-1) == Y(j-1))
            L(i,j) = L(i-1,j-1) + 1;
            b(i,j) = 3;%%%Up and left
        else
            L(i,j) = L(i-1,j-1);
        end
        if(L(i-1,j) >= L(i,j))
            L(i,j) = L(i-1,j);
            b(i,j) = 1;%Up
        end
        if(L(i,j-1) >= L(i,j))
            L(i,j) = L(i,j-1);
            b(i,j) = 2;%Left
        end
    end
end
L(:,1) = [];
L(1,:) = [];
b(:,1) = [];
b(1,:) = [];
dist = L(n,m);

D = (dist / min(m,n));
if(dist == 0)
    aLongestString = '';
else
    i = n;
    j = m;
    p = dist;
    aLongestString = {};
    while(i>0 && j>0)
        if(b(i,j) == 3)
            aLongestString{p} = X(i);
            p = p-1;
            i = i-1;
            j = j-1;
        elseif(b(i,j) == 1)
            i = i-1;
        elseif(b(i,j) == 2)
            j = j-1;
        end
    end

    if ischar(aLongestString{1})
        aLongestString = char(aLongestString)';
    else
        aLongestString = cell2mat(aLongestString);
    end
end

end