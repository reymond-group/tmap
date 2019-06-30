md build
cd build

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=%PREFIX% -DBUILD_SHARED_LIBS=ON

REM Make sure that msbuild is in the PATH and its version matches the compiler's
msbuild ALL_BUILD.vcxproj /p:Configuration=Release
msbuild INSTALL.vcxproj /p:Configuration=Release