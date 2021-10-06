@REM "%PYTHON%" setup.py install --prefix=%PREFIX%
"%PYTHON%" setup.py install
if errorlevel 1 exit 1