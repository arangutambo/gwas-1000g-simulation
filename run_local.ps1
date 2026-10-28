# Windows entry point. Runs the same bash pipeline through WSL if present, else Git Bash.
# Usage:  powershell -ExecutionPolicy Bypass -File run_local.ps1 [all|qc|tasks FROM TO|agg]
param([Parameter(ValueFromRemainingArguments = $true)] [string[]] $Args)
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $here
$argline = ($Args -join ' ')
if (Get-Command wsl -ErrorAction SilentlyContinue) {
  wsl bash -lc "cd '$(wsl wslpath -a $here)' && bash run_local.sh $argline"
} elseif (Test-Path "$env:ProgramFiles\Git\bin\bash.exe") {
  & "$env:ProgramFiles\Git\bin\bash.exe" -lc "bash run_local.sh $argline"
} else {
  Write-Error "Need WSL or Git for Windows (Git Bash). Install one, then re-run."; exit 1
}
