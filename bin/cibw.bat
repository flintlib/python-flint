rem
rem This bat file can be used to test cibuildwheel locally on Windows. The
rem cibw_*_windows.sh files have lines to set the PATH for working locally but
rem those are commented out because they are not needed in CI. To use this
rem script
rem
rem	1. Uncomment those lines
rem	2. > pip install cibuildwheel
rem	3. > bin\cibw.bat
rem
rem The variables defined below should match those that are set in CI except
rem that C:\msys64\usr\bin\bash should just be msys2 -c in CI.
rem
rem It is also worth commenting out the line to build GMP etc after you have
rem built those once because that is by far the slowest step.
rem

rem
rem If this script is run repeatedly then it would fail because of any leftover
rem wheels from a previous run so we delete them here.
rem
del /q wheelhouse\*

rem override setting in pyproject.toml
set CIBW_BUILD=cp39-* cp310-* cp311-*

cibuildwheel --platform windows
