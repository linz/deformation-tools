@echo off

for %%I in (*.in) do ..\test_okada %%I %%~dpI/out/%%~nI.out

