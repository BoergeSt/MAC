function [ M,F ] = Restrict_To_Domain( M,F,domain )

I = speye(size(M));
M(~domain,:) = I(~domain,:);
M(:,~domain) = I(:,~domain);
F(~domain) = 0;

end

