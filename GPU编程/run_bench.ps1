param(
    [int]$MinLog = 10,
    [int]$MaxLog = 22,
    [int]$Repeat = 20,
    [int]$Warmup = 3
)

$exe = Join-Path $PSScriptRoot "bin/gpu_ntt.exe"
if (-not (Test-Path $exe)) {
    Write-Host "未找到 $exe，先执行 build.bat。"
    exit 1
}

& $exe --min-log $MinLog --max-log $MaxLog --repeat $Repeat --warmup $Warmup
