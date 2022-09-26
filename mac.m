function [Ckrs] = mac(Akr,Bks)
%	function [Ckrs] = mac(Akr,Bks)
%
%	Function to compute MAC of two vector set Akr and Bks
%	
%	Example: Ckrs = mac(Akr,Bks)
%
%	by Chunyan Ao
%	2022.6.20
%   Please cite if used

[ka,r] = size(Akr);
[kb,s] = size(Bks);

if ka ~= kb
  error(['Vector size of matrix Aks (k=',num2str(ka), ... 
         ') not equal to Bkr (k=',num2str(kb),')'])
end

k = ka;

for j = 1:r
  sr(j) = (Akr(:,j)'*Akr(:,j)).^(-0.5);
end

Sr = diag(sr);

for j = 1:s
  ss(j) = (Bks(:,j)'*Bks(:,j)).^(-0.5);
end

Ss = diag(ss);

Akr = Akr * Sr;
Bks = Bks * Ss;


Dkrs = Akr' * Bks;
Ckrs = Dkrs .* Dkrs;
  