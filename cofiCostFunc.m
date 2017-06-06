function [J, grad] = cofiCostFunc(params, Y, R, num_users, num_movies, ...
                                  num_features, lambda)
%COFICOSTFUNC Collaborative filtering cost function
%   [J, grad] = COFICOSTFUNC(params, Y, R, num_users, num_movies, ...
%   num_features, lambda) returns the cost and gradient for the
%   collaborative filtering problem.
%

% Unfold the U and W matrices from params
X = reshape(params(1:num_movies*num_features), num_movies, num_features);
Theta = reshape(params(num_movies*num_features+1:end), ...
                num_users, num_features);

            
% You need to return the following values correctly
J = 0;
b1=0;
b2=0;
b3=0;
X_grad = zeros(size(X));
Theta_grad = zeros(size(Theta));
for i=1:num_movies
  for j=1:num_users
    if R(i,j)==1
      b1=(Theta(j,:)*X(i,:)'-Y(i,j)).^2;
      J=(J+b1);
          else
    end
  end
end
J=J/2;

for j=1:num_users
  for k=1:num_features
    J=J+lambda*0.5*Theta(j,k).^2;
  end
end

for i=1:num_movies
  for k=1:num_features
    J=J+lambda*0.5*X(i,k).^2;
  end
end

%for i=1:num_movies
 %for j=1:num_users
    
  %    if R(i,j)==1
   %     for k=1:num_features
   %idx=find(R(i,:)==1);
   %for l=1:idx
    %          b2=(Theta(l,:)*X(i,:)'-Y(i,l)).*Theta(l,k);
     %         X_grad(i,k)=X_grad(i,k)+b2;
  %end
   %idy=find(R(:,j)==1);
   
   %for p=1:idy  
    %         b3=(Theta(j,:)*X(p,:)'-Y(p,j)).*X(p,k);
     %        Theta_grad(j,k)=Theta_grad(j,k)+b3;
      %     end
       % end
     %end
   %end
 %end
 
% for k=1:num_features
 %for i=1:num_movies
  %idx=find(R(i,:)==1);
  %for j=1:size(idx)
  %X_grad(i,k)=sum(Theta(idx(j)).*X(i,:).*sum(Theta(idx(j),k))-sum(Y(i,idx(j)))*(sum(Theta(idx(j),k)));
 %end
%end
%end
 
 
 
 
 
 
      
      %idx=find(R(i,:)==1);
      %Thetatemp=Theta(idx,:);
      %Ytemp=Y(i,idx);
      %X_grad(i,:)=(X(i,:)*Thetatemp'-Ytemp)*Thetatemp;
      
      %idy=find(R(:,j)==1);
      %Xtemp=X(idy,:);
      %Ytemp=Y(idy,j);
      %Theta_grad(j,:)=(Xtemp'*Theta(j,:)'-Ytemp)*Xtemp;
    %end
  %end
  
  
 for k=1:num_features 
 for i=1:num_movies
  for j=1:num_users 
  if R(i,j)==1
 
     b2=b2+(Theta(j,:)*X(i,:)'-Y(i,j))*Theta(j,k);
    
  end
    
 end
 X_grad(i,k)=b2+lambda*X(i,k);
 b2=0;
 end
 end
  
  
  for k=1:num_features 
 for j=1:num_users
for i=1:num_movies
  
if R(i,j)==1

   b3=b3+(Theta(j,:)*X(i,:)'-Y(i,j))*X(i,k);
  
end
    

end
Theta_grad(j,k)=b3+lambda*Theta(j,k);
b3=0;

end
end
  
  
 
  
    
    
% for j=1:num_users 
 %for i=1:num_movies
  
  %if R(i,j)==1
    %
   %   b3=b3+(Theta(j,:)*X(i,:)'-Y(i,j));
    %end
    %for k=1:num_features
     %     Theta_grad(j,k)=b3*X(i,k);
      %  end
       % b3=0;
      %end
    %end
          
    

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost function and gradient for collaborative
%               filtering. Concretely, you should first implement the cost
%               function (without regularization) and make sure it is
%               matches our costs. After that, you should implement the 
%               gradient and use the checkCostFunction routine to check
%               that the gradient is correct. Finally, you should implement
%               regularization.
%
% Notes: X - num_movies  x num_features matrix of movie features
%        Theta - num_users  x num_features matrix of user features
%        Y - num_movies x num_users matrix of user ratings of movies
%        R - num_movies x num_users matrix, where R(i, j) = 1 if the 
%            i-th movie was rated by the j-th user
%
% You should set the following variables correctly:
%
%        X_grad - num_movies x num_features matrix, containing the 
%                 partial derivatives w.r.t. to each element of X
%        Theta_grad - num_users x num_features matrix, containing the 
%                     partial derivatives w.r.t. to each element of Theta
%
















% =============================================================

grad = [X_grad(:); Theta_grad(:)];

end
