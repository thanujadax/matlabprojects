function [VI_value,NVI] = vi(class_ind_a,class_ind_b)
%first_vector,second_cluster as [1 1 1 1 2 2 2 2 3 3 3 3 .....]
% 
% [L1m,L1n] = size(class_ind_a);
% [L2m,L2n] = size(class_ind_b);
% 
% numL1 = L1m*L1n;
% numL2 = L2m*L2n;
% 
% class_ind_a = reshape(class_ind_a,1,numL1);
% class_ind_b = reshape(class_ind_b,1,numL2);

%IMPLEMENTATION OF VARIATION OF INFORMATION
%"Comparing Clusterings",Marina Meila,2002

%Stavros Dimitriadis, 2/2008
%http://users.auth.gr/~stdimitr/index.html


%first vector
[x1 n1]=size(class_ind_a);

%second vector
[x2 n2]=size(class_ind_b);

if (n1~=n2)
    disp('two vectors dont have the same size')
end


%calculation of entropy for each cluster
%first clustering

%counting the nk (the number of point in kth cluster)

table1=tabulate(class_ind_a);
freq1=table1(:,2);
ka=length(table1(:,1));


%calculation of the entropy

entropy1=0;
for i=1:ka
    entropy1=entropy1+((freq1(i)/n1)*log10(freq1(i)/n1));
end

%second clustering
%counting the nk (the number of point in kth cluster)

table2=tabulate(class_ind_b);
freq2=table2(:,2);
kb=length(table2(:,1));


%calculation of the entropy
entropy2=0;

for i=1:kb
   entropy2=entropy2+((freq2(i)/n1)*log10(freq2(i)/n1));
end


%calling [NMI,MI,ka,kb]=MI_for_partitions(class_ind_a,class_ind_b)
[NMI,MI,ka,kb]=MI_for_partitions(class_ind_a,class_ind_b);


%define Variation of Information
entropy1 = -entropy1;
entropy2 = -entropy2;
VI_value = entropy1 + entropy2 - 2*MI ;


%normalized VI according to the paper
NVI = VI_value / log(n1);

% 
% %checking the VALUE OF mutual information according to the paper
% 
% if (MI <= min(entropy1,entropy2))
%     disp('right value of mutual information');
% end
% 
% %if two clusterings have K clusters each one ,then VI cannot be higher than
% % 2*log(K)
% no_clusters1=table1(end,1);
% no_clusters2=table2(end,1);
% 
% if (no_clusters1==no_clusters2 &&  NVI <= 2*log10(no_clusters1))
%     disp('right value of variation of information');
% end













