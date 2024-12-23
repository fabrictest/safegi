mkdir ..\obj
mkdir ..\obj\msvc10

cd ..\obj\msvc10
cmake ..\..\src -G "Visual Studio 10"
cd ..\..\src

PAUSE
