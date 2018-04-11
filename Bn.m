function [ bn ] = Bn( n, mu, mu1, N, N1, r, lambda)
% n - ������� ������������ ���������
% mu - ��������� ������������� �����
% mu1 - ��������� ������������� ����
% N - ���������� ����������� �����
% N1 - ���������� ����������� ����
% r - ������ �����
% lambfa - ����� ����� ���������

m = N1./N; % ������������� ���������� �����������
x = 2.*pi.*N.*r./lambda; % �������� ���������

b1 = mu1 .* sphbesu(n,m.*x,'j') .* (sphbesu(n,x,'j') + x .* dsphbesu(n,x,'j'));
b2 = mu .* sphbesu(n,x,'j') .* (sphbesu(n,m.*x,'j') + x.*m.*dsphbesu(n,m.*x,'j'));
b3 = mu1 .* sphbesu(n,m.*x,'j') .* (sphbesu(n,x,'h') + x .* dsphbesu(n,x,'h'));
b4 = mu .* sphbesu(n,x,'h') .* (sphbesu(n,m.*x,'j') + x.*m.*dsphbesu(n,m.*x,'j'));

bn = (b1-b2)./(b3-b4);

end

