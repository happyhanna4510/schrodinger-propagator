# =============================
# Taylor: gamma ∈ {10,20}, K ∈ {4..8}, dt ∈ {1e-4,1e-5,1e-6}
# =============================

# locate exe
$root = Split-Path -Parent $PSScriptRoot
$candidates = @(
  (Join-Path $root 'build\Release\morse.exe'),
  (Join-Path $root 'Release\morse.exe'),
  (Join-Path $root 'build\morse.exe')
)
$exe = $null; foreach ($c in $candidates) { if (Test-Path $c) { $exe = $c; break } }
if (-not $exe) { Write-Error "Не найден morse.exe."; exit 1 }
Write-Host "Использую бинарник: $exe" -ForegroundColor Green

$gammas = @(10, 20)
$Ks = 4..8
$dts = @('1e-4','1e-5','1e-6')

function Get-LogEvery([string]$dt) {
  switch ($dt) {
    '1e-4' { 10000 }
    '1e-5' { 50000 }
    '1e-6' { 100000 }
    default { 10000 }
  }
}

foreach ($g in $gammas) {
  foreach ($K in $Ks) {
    foreach ($dt in $dts) {
      $log = Get-LogEvery $dt
      Write-Host "=== Taylor: gamma=$g K=$K dt=$dt (log_every=$log) ===" -ForegroundColor Cyan
      & $exe --evolve taylor --gamma $g --K $K --dt $dt --tmax 10 --log-every $log 
    }
  }
}
