rem
rem This batch file serves the purpose of taking Windows style path arguments,
rem converting them to environment variables and calling msys2. This is needed
rem because otherwise in CI msys2 -c will mangle the paths turning e.g. C:\a\b
rem into C:ab.
rem

set tempfile=tmpfile.deleteme
set WHEELHOUSE=%1
set WHEELNAME=%2

msys2 -c bin/cibw_repair_wheel_command_windows.sh
rem C:\msys64\usr\bin\bash bin/cibw_repair_wheel_command_windows.sh
