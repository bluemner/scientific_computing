@echo off
echo ==============================================================================
echo MIT License
echo(
echo Copyright (c) 2017 Brandon Bluemner
echo( 
echo Permission is hereby granted, free of charge, to any person obtaining a copy
echo of this software and associated documentation files (the "Software"), to deal
echo in the Software without restriction, including without limitation the rights
echo to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
echo copies of the Software, and to permit persons to whom the Software is
echo furnished to do so, subject to the following conditions:
echo( 
echo The above copyright notice and this permission notice shall be included in all
echo copies or substantial portions of the Software.
echo( 
echo THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
echo IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
echo FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
echo AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
echo LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
echo OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
echo SOFTWARE.

echo ==============================================================================

set vs2013="C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\vcvarsall.bat"
if exist %vs2013% call %vs2013%  x64
if exist %vs2013% echo Visual Studio 2013 Loaded

set vs2015="C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat"
if exist %vs2015% call %vs2015% x64
if exist %vs2015% echo Visual Studio 2015 Loaded

set vs2017="C:\Program Files (x86)\Microsoft Visual Studio 16.0\VC\vcvarsall.bat"
if exist %vs2017% call %vs2017% x64
if exist %vs2017% echo Visual Studio 2017 Loaded

WHERE cl
IF %ERRORLEVEL% NEQ 0 ECHO Microsoft Visual Studio 2013 2015 not found. Visual Studio Must be 2013 or grater
IF %ERRORLEVEL% NEQ 0 Exit -1

if not exist .\bin mkdir .\bin

echo building driver
set file_name=driver
set compilerflags=/Fo.\bin\ /Od /Zi /EHsc 
set linkerflags=/OUT:bin\%file_name%.exe
cl.exe %compilerflags% source/%file_name%.cpp /link %linkerflags%

echo building house_holder_driver
set file_name=house_holder_driver
set compilerflags=/Fo.\bin\ /Od /Zi /EHsc 
set linkerflags=/OUT:bin\%file_name%.exe
cl.exe %compilerflags% source/%file_name%.cpp /link %linkerflags%
