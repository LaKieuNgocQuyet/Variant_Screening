@echo off
set "PROJECT_DIR=%~dp0"
set "VENV_DIR=%PROJECT_DIR%venv"
echo.
echo ============================================================
echo          VARIANT SCREENING PIPELINE SETUP STARTING
echo ============================================================
echo.
echo [INFO] Initializing virtual environment and installing dependencies...
echo.

python -m venv "%PROJECT_DIR%venv"
"%VENV_DIR%\Scripts\python.exe" -m pip install --upgrade pip
"%VENV_DIR%\Scripts\pip.exe" install -r "%PROJECT_DIR%src\requirements.txt"
echo.
echo ============================================================
echo        VARIANT SCREENING PIPELINE SETUP COMPLETE
echo ============================================================
echo.

pause