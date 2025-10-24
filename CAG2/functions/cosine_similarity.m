function S = cosine_similarity(A)
    epsilon = 1e-8;
    % ����ÿ�е�L2����
    norm_rows = sqrt(sum(A.^2, 2))+epsilon;
    % �����㷶�����������������㣨��ѡ��
    norm_rows(norm_rows == 0) = 1; % ��ѡ���裬������������
    % ����������
    dot_products = A * A';
    % ���㷶���˻�����
    norm_matrix = norm_rows * norm_rows';
    % �����������ƶȾ���
    S = dot_products ./ norm_matrix;
end