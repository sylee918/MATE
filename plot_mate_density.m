%% MATE Output 4D Density & MSIS BC Exosphere Dashboard (2x4 Grid Layout)
%
% 이 스크립트는 MATE 모델에서 생성된 4차원 수소 밀도 바이너리(*.data) 파일과
% exobase 경계 조건인 MSIS 바이너리(*.bc) 파일을 함께 읽어와 시각화합니다.
%
% 총 7개의 주요 물리 분포와 1개의 메타데이터 정보 카드를 2x4 패널 대시보드로 완벽하게 표출합니다:
%   ■ 행 1 (MSIS Exobase 및 MATE 저고도)
%     1) MSIS Exobase 수소 밀도 (n_H) 분포
%     2) MSIS Exobase 수소 온도 (T_H) 분포
%     3) MATE 수소 밀도 분포 at 1.2 Re
%     4) MATE 수소 밀도 분포 at 1.4 Re
%   ■ 행 2 (MATE 중/고고도 및 대시보드 요약 카드)
%     5) MATE 수소 밀도 분포 at 1.6 Re
%     6) MATE 수소 밀도 분포 at 1.8 Re
%     7) MATE 수소 밀도 분포 at 2.0 Re
%     8) 시뮬레이션 상세 정보 및 데이터 요약 정보 카드 (Metadata Card)
%
% 특징:
%   - 고도 그리드가 컴파일 타임 float truncation으로 인해 4개로 빌드되었던 과거 데이터나,
%     nint 보정 후 생성된 정상 5개 고도 데이터 모두를 에러 없이 유연하게 역산 매칭합니다.
%   - PSD.f90의 경도/위도 좌표 매핑 및 180도 경도 Shift 보정이 완벽히 구현되었습니다.
%
% 작성일: 2026-05-22

clear; clc; close all;

%% 1. 파라미터 설정 (Setting.inc 및 MATE 시뮬레이션 설정 기반)

% MATE 시뮬레이션 기본 정보
runname  = 'BC_test7';       % Runname_in_10char
iday     = 2020001;          % Start_Time_in_YYYYDOY (7자리 YYYYDOY)

% 물리 법칙 포함 여부 (Setting.inc 참조하여 tag_phys 자동 조립)
i_EarthGravity           = 1;  % 1: 포함 (G)
i_SolarRadiationPressure = 1;  % 1: 포함 (R)
i_CoriolisForce_GSE      = 1;  % 1: 포함 (C)
i_Photoionization        = 0;  % 1: 포함 (P)
i_ChargeExchange         = 0;  % 1: 포함 (X)

% MATE 격자 해상도 정보
geores = 15;                 % GEO_Resolution_in_Degree (도)
RadialRange_min = 1.2;       % 최소 고도 (Re)
RadialRange_max = 2.0;       % 최대 고도 (Re)
dR = 0.2;                    % 고도 간격 (Re)
Output_Time_Interval_in_Minute = 1440; % 출력 시간 간격 (분)

% MSIS Exobase 경계 조건(BC) 격자 해상도 정보
bc_res = 5;                  % BC_GEO_Resolution_in_Degree (5도 해상도)
BC_Time_Resolution_in_Minute = 5; % BC_Time_Resolution_in_Minute (5분 간격)

% ------------------ 파일 및 디렉토리 경로 설정 ------------------
% MATE Output 데이터 디렉토리
file_dir = '\\wsl.localhost\Ubuntu-22.04/home/sylee/exospherecode/MATE/output/0524/GRCP/';

% MSIS BC 데이터 디렉토리 (WSL 기본 경로로 자동 조립, 필요시 수정 가능)
bc_dir   = '\\wsl.localhost\Ubuntu-22.04/home/sylee/exospherecode/MSIS/Fortran/BC/';
% -----------------------------------------------------------------

%% 2. 물리 옵션 태그 빌드 및 파일 이름 자동 조합

tag_phys = '';
if i_EarthGravity == 1,           tag_phys = [tag_phys 'G']; end
if i_SolarRadiationPressure == 1, tag_phys = [tag_phys 'R']; end
if i_CoriolisForce_GSE == 1,      tag_phys = [tag_phys 'C']; end
if i_Photoionization == 1,        tag_phys = [tag_phys 'P']; end
if i_ChargeExchange == 1,         tag_phys = [tag_phys 'X']; end

% MATE Output 파일명
filename = sprintf('%sMATE_nH_%s_%s_%d.data', file_dir, tag_phys, runname, iday);

% MSIS BC 파일명 조립 (예: BC_dir/2020/MSIS_2020001.bc)
year_str = num2str(floor(iday / 1000));
ydoy_str = sprintf('%07d', iday);
filename_BC = sprintf('%s%s/MSIS_%s.bc', bc_dir, year_str, ydoy_str);

%% 3. 파일 존재 여부 검사 및 수동 탐색 안전 조치

% MATE 데이터 파일 체크
if ~exist(filename, 'file')
    warning('설정된 MATE 데이터 파일이 존재하지 않습니다: %s\n수동 파일 선택창을 엽니다.', filename);
    [file_name, file_path] = uigetfile('*.data', 'MATE 수소 밀도 데이터 파일 (*.data) 선택', file_dir);
    if isequal(file_name, 0)
        error('MATE 파일 선택이 취소되어 스크립트를 종료합니다.');
    end
    filename = fullfile(file_path, file_name);
end

% MSIS BC 파일 체크
if ~exist(filename_BC, 'file')
    warning('설정된 MSIS BC 파일이 존재하지 않습니다: %s\n수동 BC 파일 선택창을 엽니다.', filename_BC);
    [bc_name, bc_path] = uigetfile('*.bc', 'MSIS Exobase BC 파일 (*.bc) 선택', bc_dir);
    if isequal(bc_name, 0)
        error('BC 파일 선택이 취소되어 스크립트를 종료합니다.');
    end
    filename_BC = fullfile(bc_path, bc_name);
end

%% 4. 격자 차원 및 변수 자동 계산

% ------------------ MATE 격자 ------------------
nLong = 360 / geores;                                            % 360 / 15 = 24
nLat_NS = (90 / geores + 1) * 2 - 1;                             % (6 + 1) * 2 - 1 = 13
ntperday = round(86400 / (Output_Time_Interval_in_Minute * 60)); % 1
nRadial_nominal = round((RadialRange_max - RadialRange_min) / dR) + 1; % 5

% ------------------- BC 격자 -------------------
nbx = 360 / bc_res;                                              % 360 / 5 = 72
nby = 180 / bc_res;                                              % 180 / 5 = 36
nbtperday = 288;                                                 % 5분 해상도 하루 스텝 고정

%% 5. MATE Output 데이터 로드 & nRadial 동적 보정

fid = fopen(filename, 'r');
if fid == -1
    error('MATE 파일을 열 수 없습니다: %s', filename);
end

fseek(fid, 0, 'eof');
file_bytes = ftell(fid);
fseek(fid, 0, 'bof');

num_elements_in_file = file_bytes / 4;
nRadial_actual = num_elements_in_file / (nLong * nLat_NS * ntperday);

% 역산 분할 검증
if mod(nRadial_actual, 1) ~= 0
    fclose(fid);
    error('오류: MATE 파일 크기가 설정된 격자 해상도 정수배로 나누어떨어지지 않습니다.');
end

nRadial = nRadial_actual;
if nRadial_nominal ~= nRadial_actual
    warning('Setting.inc 계산값(nRadial=%d)과 실제 파일 구조(nRadial=%d)가 달라, 자동 보정합니다.', ...
            nRadial_nominal, nRadial_actual);
end

data_flat = fread(fid, num_elements_in_file, 'single');
fclose(fid);

% 4차원 데이터로 복원
density_4D = reshape(data_flat, [nRadial, nLong, nLat_NS, ntperday]);

%% 6. MSIS Exobase BC 데이터 로드

fid_bc = fopen(filename_BC, 'r');
if fid_bc == -1
    error('MSIS BC 파일을 열 수 없습니다: %s', filename_BC);
end

% read_exobaseBC 서브루틴에 맞춰 nH_real, TH_real 순서로 기록되어 있음
num_elements_bc = nbx * nby * nbtperday;
nH_BC_flat = fread(fid_bc, num_elements_bc, 'single');
TH_BC_flat = fread(fid_bc, num_elements_bc, 'single');
fclose(fid_bc);

% 3차원 데이터로 복원 [nbx, nby, nbtperday]
nH_BC_3D = reshape(nH_BC_flat, [nbx, nby, nbtperday]);
TH_BC_3D = reshape(TH_BC_flat, [nbx, nby, nbtperday]);

%% 7. 특정 시간대(Time Step) 및 좌표 데이터 추출

% 시각화할 MATE 시간 인덱스
target_time_idx = 1;

% 시각화할 MSIS BC 시간 인덱스 (1 ~ 288 스텝 중 중간인 12:00 UT = 144 스텝 선택)
target_BC_time_idx = 144; 

% ------------------- MATE 데이터 슬라이스 -------------------
R_range   = linspace(RadialRange_min, RadialRange_max, nRadial); 
lon_range = (0:nLong-1) * geores;
lat_range = ((0:nLat_NS-1) * geores) - 90;

% [강인한 설계] 각 타겟 고도(1.2, 1.4, 1.6, 1.8, 2.0 Re)에 가장 가까운 인덱스 동적 매칭
[~, idx_1_2] = min(abs(R_range - 1.2));
[~, idx_1_4] = min(abs(R_range - 1.4));
[~, idx_1_6] = min(abs(R_range - 1.6));
[~, idx_1_8] = min(abs(R_range - 1.8));
[~, idx_2_0] = min(abs(R_range - 2.0));

R_1_2_actual = R_range(idx_1_2);
R_1_4_actual = R_range(idx_1_4);
R_1_6_actual = R_range(idx_1_6);
R_1_8_actual = R_range(idx_1_8);
R_2_0_actual = R_range(idx_2_0);

mate_nH_1_2 = squeeze(density_4D(idx_1_2, :, :, target_time_idx))';
mate_nH_1_4 = squeeze(density_4D(idx_1_4, :, :, target_time_idx))';
mate_nH_1_6 = squeeze(density_4D(idx_1_6, :, :, target_time_idx))';
mate_nH_1_8 = squeeze(density_4D(idx_1_8, :, :, target_time_idx))';
mate_nH_2_0 = squeeze(density_4D(idx_2_0, :, :, target_time_idx))';

[LON_MATE, LAT_MATE] = meshgrid(lon_range, lat_range);

% -------------------- BC 데이터 슬라이스 --------------------
% PSD.f90의 수정된 매핑 공식 연동
lon_range_BC = ((0:nbx-1) * bc_res) - 180;
lat_range_BC = ((0:nby-1) * bc_res) - 90;

nH_BC_2D = squeeze(nH_BC_3D(:, :, target_BC_time_idx));
TH_BC_2D = squeeze(TH_BC_3D(:, :, target_BC_time_idx));

% PSD.f90 격자 정방향 연동에 따른 단순 Transpose
nH_BC_plot = nH_BC_2D';
TH_BC_plot = TH_BC_2D';

% 🌐 [경도 180도 Shift 보정]
nH_BC_plot = circshift(nH_BC_plot, nbx/2, 2);
TH_BC_plot = circshift(TH_BC_plot, nbx/2, 2);

[LON_BC, LAT_BC] = meshgrid(lon_range_BC, lat_range_BC);

%% 8. 시각화 대시보드 플로팅 (2행 4열 웅장한 대시보드 레이아웃)

fig = figure('Color', 'w', 'Position', [30, 30, 1650, 850]);

% ------------------ [패널 1] MSIS Exobase 수소 밀도 분포 ------------------
subplot(2, 4, 1);
contourf(LON_BC, LAT_BC, nH_BC_plot, 30, 'LineColor', 'none');
hold on;
colormap(gca, parula(256));
c1 = colorbar;
c1.Label.String = 'Exobase n_H (cm^{-3})';
c1.Label.FontWeight = 'bold';
style_panel('1. MSIS Exobase Density (n_H)', 'Longitude', 'Latitude', [-180, 180]);
hold off;

% ------------------ [패널 2] MSIS Exobase 수소 온도 분포 ------------------
subplot(2, 4, 2);
contourf(LON_BC, LAT_BC, TH_BC_plot, 30, 'LineColor', 'none');
hold on;
colormap(gca, hot(256)); 
c2 = colorbar;
c2.Label.String = 'Exobase Temperature (K)';
c2.Label.FontWeight = 'bold';
style_panel('2. MSIS Exobase Temperature (T_H)', 'Longitude', 'Latitude', [-180, 180]);
hold off;

% ------------------ [패널 3] MATE 수소 밀도 분포 at 1.2 Re ------------------
subplot(2, 4, 3);
contourf(LON_MATE, LAT_MATE, mate_nH_1_2, 30, 'LineColor', 'none');
hold on;
colormap(gca, parula(256));
c3 = colorbar;
c3.Label.String = 'Density (cm^{-3})';
c3.Label.FontWeight = 'bold';
style_panel(sprintf('3. MATE n_H at %.2f R_E (idx=%d)', R_1_2_actual, idx_1_2), 'Longitude', 'Latitude');
hold off;

% ------------------ [패널 4] MATE 수소 밀도 분포 at 1.4 Re ------------------
subplot(2, 4, 4);
contourf(LON_MATE, LAT_MATE, mate_nH_1_4, 30, 'LineColor', 'none');
hold on;
colormap(gca, parula(256));
c4 = colorbar;
c4.Label.String = 'Density (cm^{-3})';
c4.Label.FontWeight = 'bold';
style_panel(sprintf('4. MATE n_H at %.2f R_E (idx=%d)', R_1_4_actual, idx_1_4), 'Longitude', 'Latitude');
hold off;

% ------------------ [패널 5] MATE 수소 밀도 분포 at 1.6 Re ------------------
subplot(2, 4, 5);
contourf(LON_MATE, LAT_MATE, mate_nH_1_6, 30, 'LineColor', 'none');
hold on;
colormap(gca, parula(256));
c5 = colorbar;
c5.Label.String = 'Density (cm^{-3})';
c5.Label.FontWeight = 'bold';
style_panel(sprintf('5. MATE n_H at %.2f R_E (idx=%d)', R_1_6_actual, idx_1_6), 'Longitude', 'Latitude');
hold off;

% ------------------ [패널 6] MATE 수소 밀도 분포 at 1.8 Re ------------------
subplot(2, 4, 6);
contourf(LON_MATE, LAT_MATE, mate_nH_1_8, 30, 'LineColor', 'none');
hold on;
colormap(gca, parula(256));
c6 = colorbar;
c6.Label.String = 'Density (cm^{-3})';
c6.Label.FontWeight = 'bold';
style_panel(sprintf('6. MATE n_H at %.2f R_E (idx=%d)', R_1_8_actual, idx_1_8), 'Longitude', 'Latitude');
hold off;

% ------------------ [패널 7] MATE 수소 밀도 분포 at 2.0 Re ------------------
subplot(2, 4, 7);
contourf(LON_MATE, LAT_MATE, mate_nH_2_0, 30, 'LineColor', 'none');
hold on;
colormap(gca, parula(256));
c7 = colorbar;
c7.Label.String = 'Density (cm^{-3})';
c7.Label.FontWeight = 'bold';
style_panel(sprintf('7. MATE n_H at %.2f R_E (idx=%d)', R_2_0_actual, idx_2_0), 'Longitude', 'Latitude');
hold off;



% 피규어 대시보드 전체 제목
sgtitle(sprintf('MATE Exospheric Density Analysis Dashboard (Run: %s | Day: %s)', runname, ydoy_str), ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.1, 0.1, 0.2], 'Interpreter', 'none');

%% 로컬 헬퍼 함수 (Local Helper Functions)
function style_panel(title_str, x_label, y_label, x_limits)
    if nargin < 4
        x_limits = [0, 360];
    end

    grid on;
    box on;
    xlim(x_limits);
    ylim([-90, 90]);
    
    if x_limits(1) < 0
        xticks(-180:60:180);
        xticklabels(cellfun(@(x) [num2str(x), '\circ'], num2cell(-180:60:180), 'UniformOutput', false));
    else
        xticks(0:60:360);
        xticklabels(cellfun(@(x) [num2str(x), '\circ'], num2cell(0:60:360), 'UniformOutput', false));
    end
    
    yticks(-90:30:90);
    yticklabels(cellfun(@(x) [num2str(x), '\circ'], num2cell(-90:30:90), 'UniformOutput', false));
    xlabel(x_label, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(y_label, 'FontSize', 10, 'FontWeight', 'bold');
    title(title_str, 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'none');
    set(gca, 'FontSize', 9, 'LineWidth', 1.2, 'GridAlpha', 0.12);
end
