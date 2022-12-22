SET
FOR /F "tokens=* USEBACKQ" %%F IN (`python -c "import sys, os; print(os.path.dirname(sys.executable))"`) DO (
	SET VIRTUAL_ENV_BIN=%%F
)
echo %VIRTUAL_ENV_BIN%
msys2 -c bin/cibw_before_build_windows.sh
