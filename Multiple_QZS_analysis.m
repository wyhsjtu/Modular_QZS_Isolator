clear all; clc;
color_full = ["#5f5f5f",  "#7262ac","#2e7ebb" ,"#2e974e" ,"#e25508","#d92523"   ];
colors_2 = ["#cecece",  "#cfcfe5","#b7d4ea","#b8e3b2","#fdc38d" ,"#fcab8f"]; 

%% 1014- 1015 cases
%--------------- Case 1: 70mm ---------------
amp_70 = 1;

filename_70mm = 'Multiple_QZS/Config2_2.xlsx';
data_70mm = readtable(filename_70mm, 'Sheet', 'Config2_2');
load_70mm = data_70mm.LoadValue*amp_70;
defl_70mm = data_70mm.PositionValue;

%--------------- Case 2: 50mm ---------------
amp_50 = 1;

filename_50mm = 'Multiple_QZS/Config3_1.xlsx';
data_50mm = readtable(filename_50mm, 'Sheet', 'Config3_1');
load_50mm = data_50mm.LoadValue*amp_50;
defl_50mm = data_50mm.PositionValue;

%--------------- Case 3: 60mm ---------------
amp_60 = 1.6;

filename_60mm = 'Multiple_QZS/Config1_2.xlsx';
data_60mm = readtable(filename_60mm, 'Sheet', 'Config1_2');
load_60mm = data_60mm.LoadValue*amp_60;
defl_60mm = data_60mm.PositionValue*1.15;

%--------------- Plot ---------------
figure()
grid on
hold on
plot(defl_50mm,load_50mm,'Color', color_full(5), 'LineWidth', 1.2);
plot(defl_60mm,load_60mm,'Color', color_full(3), 'LineWidth', 1.2);
plot(defl_70mm,load_70mm,'Color', color_full(4), 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('50mm','60mm','70mm','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('Multiple Configuration QZS', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')

hold off
box on
%% 1017 cases
%--------------- Case 1: 50mm ---------------
amp_50 = 1;

file_50mm = 'Multiple_QZS/1017_50mm_1.xlsx';
data_n_50mm = readtable(file_50mm, 'Sheet', '1017_50mm_1');
load_n_50mm = data_n_50mm.LoadValue*amp_50;
defl_n_50mm = data_n_50mm.PositionValue;

%--------------- Case 2: 70mm ---------------
amp_70 = 1;

file_n_70mm = 'Multiple_QZS/1017_70mm_1.xlsx';
data_n_70mm = readtable(file_n_70mm, 'Sheet', '1017_70mm_1');
load_n_70mm = data_n_70mm.LoadValue*amp_50;
defl_n_70mm = data_n_70mm.PositionValue;



figure()
grid on
hold on
plot(defl_n_50mm,load_n_50mm,'Color', color_full(5), 'LineWidth', 1.2);
plot(defl_n_70mm,load_n_70mm,'Color', color_full(3), 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('50mm','70mm','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('Multiple Configuration QZS, 50mm & 70mm', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
box on

%% Multiple combination

file_comb_1 = 'Multiple_QZS/5070_combination1_1.xlsx';
data_comb_1 = readtable(file_comb_1, 'Sheet', '5070_combination1_1');
load_comb_1 = data_comb_1.LoadValue;
defl_comb_1 = data_comb_1.PositionValue;

file_comb_2 = 'Multiple_QZS/5070_combination1_2.xlsx';
data_comb_2 = readtable(file_comb_2, 'Sheet', '5070_combination1_2');
load_comb_2 = data_comb_2.LoadValue;
defl_comb_2 = data_comb_2.PositionValue;

% figure()
% grid on
% hold on
% plot(defl_comb_1,load_comb_1,'Color', color_full(5), 'LineWidth', 1.2);
% plot(defl_comb_2,load_comb_2,'Color', color_full(3), 'LineWidth', 1.2);
% xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
% ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
% legend('Position 1','Position 2','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
% title('QZS Combination 1, two directions', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
% hold off
% box on

%------------------------------------------------
file_comb_3 = 'Multiple_QZS/5070_combination2_1.xlsx';
data_comb_3 = readtable(file_comb_3, 'Sheet', '5070_combination2_1');
load_comb_3 = data_comb_3.LoadValue;
defl_comb_3 = data_comb_3.PositionValue;

file_comb_4 = 'Multiple_QZS/5070_combination2_2.xlsx';
data_comb_4 = readtable(file_comb_4, 'Sheet', '5070_combination2_2');
load_comb_4 = data_comb_4.LoadValue;
defl_comb_4 = data_comb_4.PositionValue;

figure()
grid on
hold on
plot(defl_comb_3,load_comb_3,'Color', color_full(5), 'LineWidth', 1.2);
plot(defl_comb_4,load_comb_4,'Color', color_full(3), 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('Position 1','Position 2','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('QZS Combination 2, two directions', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
box on

%------------------------------------------------
file_comb_5 = 'Multiple_QZS/5070_combination3_1.xlsx';
data_comb_5 = readtable(file_comb_5, 'Sheet', '5070_combination3_1');
load_comb_5 = data_comb_5.LoadValue;
defl_comb_5 = data_comb_5.PositionValue;

file_comb_6 = 'Multiple_QZS/5070_combination3_2.xlsx';
data_comb_6 = readtable(file_comb_6, 'Sheet', '5070_combination3_2');
load_comb_6 = data_comb_6.LoadValue;
defl_comb_6 = data_comb_6.PositionValue;

figure()
grid on
hold on
plot(defl_comb_5,load_comb_5,'Color', color_full(5), 'LineWidth', 1.2);
plot(defl_comb_6,load_comb_6,'Color', color_full(3), 'LineWidth', 1.2);
xlabel('Displacement [mm]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
ylabel('Loading [N]', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
legend('Position 1','Position 2','FontName', 'Calibri', 'FontSize', 12, 'FontWeight', 'bold')
title('QZS Combination 3, two directions', 'FontName', 'Calibri', 'FontSize', 15, 'FontWeight', 'bold')
hold off
box on

