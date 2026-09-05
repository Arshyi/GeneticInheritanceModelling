param(
    [ValidateSet('test','reproduce','science','benchmark','manuscript','fetch')]
    [string]$Stage='test',
    [string]$PythonExecutable='C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    [Parameter(ValueFromRemainingArguments=$true)]
    [string[]]$StageArguments
)
$ErrorActionPreference='Stop'
Push-Location -LiteralPath $PSScriptRoot
try {
    & $PythonExecutable -X utf8 run.py $Stage @StageArguments
    if ($LASTEXITCODE -ne 0) { throw "Research stage failed with exit code $LASTEXITCODE" }
}
finally { Pop-Location }
