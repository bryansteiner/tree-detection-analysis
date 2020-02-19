X1 = transpose(xlsread('data/treeinfo.xlsx', 1, 'B3:BX29'));
Y1 = transpose(xlsread('data/treeinfo.xlsx', 1, 'B30:BX30'));

B = mnrfit(X1,Y1);
