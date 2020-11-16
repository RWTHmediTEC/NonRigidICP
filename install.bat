@echo off

echo Library installer script
echo ------------------------
echo[

set LIBRARY_NAME=NonRigidICP
set VERSION=1.8.1

if "%MATLAB_LIBS%" == "" (
    echo Error: The environment variable 'MATLAB_LIBS' is not set.
    pause
    exit -1
) else (
    echo Installing "%LIBRARY_NAME%-%VERSION%" in "%MATLAB_LIBS%".
    echo Close this window to abort the installation.
    timeout /t 15

    if not exist "%MATLAB_LIBS%" (
        echo Creating folder "%MATLAB_LIBS%"...
        mkdir "%MATLAB_LIBS%"
    )

    if not exist "%MATLAB_LIBS%\%LIBRARY_NAME%-%VERSION%" (
        echo Creating folder "%LIBRARY_NAME%-%VERSION%" in "%MATLAB_LIBS%"...
        mkdir "%MATLAB_LIBS%\%LIBRARY_NAME%-%VERSION%"
    )

    echo Copying the data...
    if exist exclude.txt rm exclude.txt
    echo %~dp0.git >> exclude.txt
    echo %~dp0src\codegen >> exclude.txt
    echo exclude.txt >> exclude.txt
    xcopy . "%MATLAB_LIBS%\%LIBRARY_NAME%-%VERSION%\" /f /y /e /EXCLUDE:exclude.txt
    rm exclude.txt

    pause
)
