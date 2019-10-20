% Task: check whether A, B is in canonical form, and return p, Index if yes

function [cano_flag, Index, p, n,m,pi] = check_cano(A, B)


% cano_flag: if 1, then is canonical, if 0, then not

cano_flag=1;


p=0;
pi=0;

% check A, B's size
n = length(A);
m = size(B,2);
if size(A) ~= [n,n] % A should be n by n
    cano_flag=0;
    
end

if size(B,1)~= n
    cano_flag =0; % B should be n by m
    
end
Index = zeros(m,1);

%% check B 
% Check B
% 1) B's column should have one and only one 1, others are 0
% 2) B's 1 location should increase
% 3) B's last column's 1 should be at n


Kk = zeros(m+1,0); % Kk = [k0,..., km], index = Kk(2:end)
for i =1:m
    % check B
    ki= find(B(:,i)~=0);
    if length(ki)~= 1  % if have more than one entries nonzero
        cano_flag =0;
        disp("B not cano: B's colum has more than one entries nonzero")
        break
    end
    if B(ki,i)~=1  % if the only one nonzero entry not 1
        cano_flag =0;
        disp("B not cano: B's colum's non zero is not 1")
        break
    end
    Kk(i+1)=ki; 
    if Kk(i+1)<= Kk(i)  % if B's ones' row index not increase
        cano_flag=0;
        disp("B not cano: B's ki index not increase")
        break
    end
    
    
end
if cano_flag ==1
    if B(end,end)~=1
        cano_flag =0;
        disp("B not cano: B's last colum's 1 is not at the end")
    end
end


%% Check A
% at B's zero rows, A(ki, ki+1)=1, others 0

if cano_flag ==1 % if B is cano
for i =1:n
    Ai = A(i,:);
    chek = find(Kk == i); % check if row i is in B's nonzero entry
    if length(chek) == 0 % if i is not B's nonzero entry
        Ai_nonzero = find(Ai ~=0);
        
        if length(Ai_nonzero)~=1 % if not only one nonzero entry
            cano_flag =0;
            disp("A not cano: length(Ai_nonzero)~=1 ")
            break
        end
        
        if Ai_nonzero ~= i+1  % if A's nonzero entry is not at i+1
            cano_flag =0;
            disp("A not cano: Ai_nonzero ~= i+1 ")
            break
        end
        
        if A(i, Ai_nonzero)~=1 % if A's nonzero entry is not 1
            cano_flag=0;
            disp("A not cano: A(i, Ai_nonzero)~=1")
            break
        end
    end
end
end

        
        
 if cano_flag==1
     Index = Kk(2:end);
     pi =Index -Kk(1:end-1);
     p = max(pi);
 end
 





