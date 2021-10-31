@REM "%PYTHON%" setup.py install --prefix=%PREFIX%
%PYTHON% -m pip install . -vv --no-build-isolation
@REM "%PYTHON%" setup.py install

@REM if errorlevel 1 exit 1