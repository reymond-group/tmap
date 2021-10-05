SET CMAKE_GENERATOR="Visual Studio 14"
"%PYTHON%" setup.py install --prefix=%PREFIX%
REM "%PYTHON%" setup.py install
if errorlevel 1 exit 1