
cl /c /I. /I"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\include" /Foutil.obj util.c
cl /c /EHsc /I. /I"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\include" /Fomain.obj main.cpp

cl /Femain.exe main.obj util.obj /link
