@echo off
title Variant Screening Launcher
chcp 65001 > nul
setlocal enabledelayedexpansion

set OUTPUT_DIR=OUTPUT
set "PROJECT_DIR=%~dp0"
set "VENV_DIR=%PROJECT_DIR%venv"

echo.
echo ============================================================
echo               VARIANT SCREENING PIPELINE LAUNCHER
echo ============================================================
echo.

set /p INPUT_DIR= Enter input directory path:
set /p GENE_NAME= Enter gene name(s), separated by comma:
set /p OUTPUT_FILE_NAME= Enter output file name (no extension):

set INPUT_DIR=%INPUT_DIR:"=%
set GENE_NAME=%GENE_NAME:"=%
set OUTPUT_FILE_NAME=%OUTPUT_FILE_NAME:"=%

if not exist "!INPUT_DIR!" (
    echo [ERROR] Input directory not found: "!INPUT_DIR!"
    pause
    exit /b
)

echo.
echo ------------------------------------------------------------

"%VENV_DIR%\Scripts\python.exe" "%~dp0src\VariantScreening_v0.2.py" -I "!INPUT_DIR!" --gene_name "!GENE_NAME!" -O "%OUTPUT_DIR%\%OUTPUT_FILE_NAME%.xlsx"

pause