function [ephysKilosortPath, ephysRawDir, ephysMetaDir, saveLocation, savePath, ...
    decompressDataLocal] = bc_definePaths

ephysKilosortPath = 'C:\Data\FG001\2022-06-01';% path to your kilosort output files 
ephysRawDir = dir('C:\Data\FG001\2022-06-01\2022-06-01 FG001_g1_t0.imec0.ap.bin'); % path to yourraw .bin or .dat data
ephysMetaDir = dir('C:\Data\FG001\2022-06-01\2022-06-01 FG001_g1_t0.imec0.ap.meta'); % path to your .meta or .oebin meta file
saveLocation = 'C:\Data\FG001\2022-06-01\bombcell'; % where you want to save the quality metrics 
savePath = fullfile(saveLocation, 'qMetrics'); 
decompressDataLocal = 'C:\Data\FG001\2022-06-01'; % where to save raw decompressed ephys data 
% 
% ephysKilosortPath = 'E:\FG005\2023-08-21\FG005_2023-08-21_g0\FG005_2023-08-21_g0_imec0\phy_curated';% path to your kilosort output files 
% ephysRawDir = dir('E:\FG005\2023-08-21\FG005_2023-08-21_g0\FG005_2023-08-21_g0_imec0\FG005_2023-08-21_g0_t0.imec0.ap.bin'); % path to yourraw .bin or .dat data
% ephysMetaDir = dir('E:\FG005\2023-08-21\FG005_2023-08-21_g0\FG005_2023-08-21_g0_imec0\FG005_2023-08-21_g0_t0.imec0.ap.meta'); % path to your .meta or .oebin meta file
% saveLocation = 'E:\FG005\2023-08-21\FG005_2023-08-21_g0\FG005_2023-08-21_g0_imec0\bombcell_new'; % where you want to save the quality metrics 
% savePath = fullfile(saveLocation, 'qMetrics'); 
% decompressDataLocal = '/media/julie/ExtraHD/decompressedData'; % where to save raw decompressed ephys data 