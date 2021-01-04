%% 导入数据
clear
clc
load('doppler_01.mat');% QR0 QR3 真值
load('SDTest_Data.mat');% 读取的多普勒单差
load('DDTest_Data.mat');% 读取的多普勒双差
% 单差
load(GetStructName('doppler_all'));% matlab仿真出的多普勒
Single_Diff(Fixed_Data, True_Data);
% 双差
load(GetStructName('RealSingleDiff_Data'));% matlab仿真出的单差
Double_Diff(RealSingleDiff_Data, 'G19', 'C30');
load(GetStructName('DoubleDiff_Data'));% matlab仿真出的双差

%% 与接收机中的多普勒对比
sate_N = length(Fixed_Data.PRN);
for i = 1:sate_N
    if length(Fixed_Data.time{i})<100
        continue
    end
    
    Prn_Index = strcmp(QR3_Data.PRN, Fixed_Data.PRN{i});
    if sum(Prn_Index) == 0
        continue;
    elseif sum(Prn_Index) == 1
        myplot1(Fixed_Data.time{i}, Fixed_Data.doppler{i}, 'doppler', Fixed_Data.PRN{i}, i, 'r');
        hold on
        myplot1(QR3_Data.time{Prn_Index}, QR3_Data.doppler{Prn_Index}, 'doppler', Fixed_Data.PRN{i}, i, 'blue');
        hold on
        legend('接收机','仿真');
    end
    
end
%% 真值
Sate_N = length(True_Data.PRN);

for i = 1:Sate_N
    myplot1(True_Data.time{i}, True_Data.doppler{i}, 'doppler', True_Data.PRN{i}, i, 'r');
end


%% plot 单差对比

plotcompare(RealSingleDiff_Data, SDTest_Data);
%% plot 双差对比
plotcompare(DoubleDiff_Data, DDTest_Data);
%% plot双差
Sate_N = length(DoubleDiff_Data.PRN);

for i = 1:Sate_N
    myplot1(DoubleDiff_Data.time{i}, DoubleDiff_Data.doppler{i}, 'doppler', DoubleDiff_Data.PRN{i}, i, 'r');
end
%% 测试
Data = RealSingleDiff_Data;
ref_gps = 'G19';
ref_bd = 'C30';
PRNtmp = {};
timetemp = {};
doptemp = {};
k = 1;
if sum(strcmp(Data.PRN, ref_gps)) == 0
    errordlg('选中的参考星无效');
    return;
elseif sum(strcmp(Data.PRN, ref_gps)) == 1
    ref_gps_index = find(strcmp(Data.PRN, ref_gps));
end
if sum(strcmp(Data.PRN, ref_bd)) == 0
    errordlg('选中的参考星无效');
    return;
elseif sum(strcmp(Data.PRN, ref_bd)) == 1
    ref_bd_index = find(strcmp(Data.PRN, ref_bd));
end

for i = 1:length(Data.PRN)
    if strcmp(Data.PRN{i}, ref_gps)
        continue;
    end
    
    if strcmp(Data.PRN{i}, ref_bd)
        continue;
    end
    doptemp{k} = [];
    timetemp{k} = [];
    if Data.PRN{i}(1) == 'G'
        PRNtmp{k} = [Data.PRN{i} '_' ref_gps(2:end)];
        for j = 1:length(Data.time{i})
            indextmp = (Data.time{ref_gps_index} == Data.time{i}(j));
            if sum(indextmp) == 0
                continue;
            elseif sum(indextmp) == 1
                doptemp{k} = [doptemp{k} Data.doppler{i}(j) - Data.doppler{ref_gps_index}(indextmp)];
                timetemp{k} = [timetemp{k} Data.time{i}(j)];
            end
        end
        k = k+1;
    elseif Data.PRN{i}(1) == 'C'            
        PRNtmp{k} = [Data.PRN{i} '_' ref_bd(2:end)];
        for j = 1:length(Data.time{i})                
            indextmp = (Data.time{ref_bd_index} == Data.time{i}(j));
            if sum(indextmp) == 0
                continue;
            elseif sum(indextmp) == 1
                doptemp{k} = [doptemp{k} Data.doppler{i}(j) - Data.doppler{ref_bd_index}(indextmp)];
                timetemp{k} = [timetemp{k} Data.time{i}(j)];
            end
        end
        k = k+1;
    end
end
DoubleDiff_Data = GetDataStruct(timetemp, doptemp, PRNtmp);
save('DoubleDiff_Data.mat', 'DoubleDiff_Data')
%% 函数定义
% 画图
function myplot1(X, Y, ylabel_s, title_str, n, color)
%PLOTWAV 此处显示有关此函数的摘要
%   此处显示详细说明

figure(n);set(gcf,'Position',get(0,'ScreenSize'));
plot(X, Y, 'Color', color);
axis xy
axis tight
ylabel(ylabel_s);
xlabel('Time (secs)');
title(title_str);
end

% 画图对比 
function plotcompare(Data1, Data2)
Sate_N = length(Data1.PRN);
for i = 1:Sate_N
    if length(Data1.time{i})<100
        continue
    end
    
    Prn_Index = strcmp(Data2.PRN, Data1.PRN{i});
    if sum(Prn_Index) == 0
        continue;
    elseif sum(Prn_Index) == 1
        myplot1(Data1.time{i}, Data1.doppler{i}, 'doppler', strrep(Data1.PRN{i}, '_', '\_'), i, 'r');
        hold on
        %myplot1(Data2.time{Prn_Index}, Data2.doppler{Prn_Index}, 'doppler', strrep(Data1.PRN{i}, '_', '\_'), i, 'blue');
        plot(Data2.time{Prn_Index}, Data2.doppler{Prn_Index}, '--', 'Color', 'blue');
        axis xy
        axis tight
        ylabel('doppler');
        xlabel('Time (secs)');
        hold off
        legend('data1','data2');
    end
end
end
% 单差
function Single_Diff(Data, Ref_Data)
    PRNtmp = {};
    timetmp = {};
    doptmp = {};

    N = length(Data.PRN);
    for i = 1:N
        index = strcmp(Ref_Data.PRN, Data.PRN{i});
        if sum(index) == 0
            continue;
        elseif sum(index) == 1
            time_N = length(Data.time{i});
            for j = 1:time_N
                time_index = (Ref_Data.time{index} == Data.time{i}(j));
                if sum(time_index) == 0
                    continue;
                elseif sum(time_index) == 1
                    SingleDiffDop = Data.doppler{i}(j) - Ref_Data.doppler{index}(time_index);
                    % 存数据
                    index_tmp = strcmp(PRNtmp, Data.PRN{i});
                    if sum(index_tmp) == 0
                        PRNtmp = [PRNtmp, Data.PRN{i}];
                        timetmp = [timetmp, Data.time{i}(j)];
                        doptmp = [doptmp, SingleDiffDop];
                    elseif sum(index_tmp) == 1
                        timetmp{index_tmp} = [timetmp{index_tmp}, Data.time{i}(j)];
                        doptmp{index_tmp} = [doptmp{index_tmp}, SingleDiffDop];
                    end
                end
            end
        end
    end
    RealSingleDiff_Data = GetDataStruct(timetmp, doptmp, PRNtmp);
    save(GetStructName('RealSingleDiff_Data'), 'RealSingleDiff_Data');
    
end
% 双差
function Double_Diff(Data, ref_gps, ref_bd)    
PRNtmp = {};
timetemp = {};
doptemp = {};
k = 1;
if sum(strcmp(Data.PRN, ref_gps)) == 0
    errordlg('选中的参考星无效');
    return;
elseif sum(strcmp(Data.PRN, ref_gps)) == 1
    ref_gps_index = find(strcmp(Data.PRN, ref_gps));
end
if sum(strcmp(Data.PRN, ref_bd)) == 0
    errordlg('选中的参考星无效');
    return;
elseif sum(strcmp(Data.PRN, ref_bd)) == 1
    ref_bd_index = find(strcmp(Data.PRN, ref_bd));
end

for i = 1:length(Data.PRN)
    if strcmp(Data.PRN{i}, ref_gps)
        continue;
    end
    
    if strcmp(Data.PRN{i}, ref_bd)
        continue;
    end
    doptemp{k} = [];
    timetemp{k} = [];
    if Data.PRN{i}(1) == 'G'
        PRNtmp{k} = [Data.PRN{i} '_' ref_gps(2:end)];
        for j = 1:length(Data.time{i})
            indextmp = (Data.time{ref_gps_index} == Data.time{i}(j));
            if sum(indextmp) == 0
                continue;
            elseif sum(indextmp) == 1
                doptemp{k} = [doptemp{k} Data.doppler{i}(j) - Data.doppler{ref_gps_index}(indextmp)];
                timetemp{k} = [timetemp{k} Data.time{i}(j)];
            end
        end
        k = k+1;
    elseif Data.PRN{i}(1) == 'C'            
        PRNtmp{k} = [Data.PRN{i} '_' ref_bd(2:end)];
        for j = 1:length(Data.time{i})                
            indextmp = (Data.time{ref_bd_index} == Data.time{i}(j));
            if sum(indextmp) == 0
                continue;
            elseif sum(indextmp) == 1
                doptemp{k} = [doptemp{k} Data.doppler{i}(j) - Data.doppler{ref_bd_index}(indextmp)];
                timetemp{k} = [timetemp{k} Data.time{i}(j)];
            end
        end
        k = k+1;
    end
end
DoubleDiff_Data = GetDataStruct(timetemp, doptemp, PRNtmp);
save(GetStructName('DoubleDiff_Data'), 'DoubleDiff_Data')
end
function [Data] = GetDataStruct(time, dop, PRN)
Data.time = time;
Data.PRN =PRN;
Data.doppler = dop;
end

function StructName = GetStructName(NameStr)
load('/home/ftw/work_space/Matlab/Gittest/Branch.mat');
StructName = [NameStr, '_', num2str(Branch), '.mat'];
end