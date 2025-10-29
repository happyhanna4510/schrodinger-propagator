# ============================================
# Taylor runs in parallel (WinPS 5.1 & PS7 OK)
# - Static Morse written once per gamma -> results/morse/g{gamma}/
# - Taylor results -> results/taylor/g{gamma}/K{K}/dt_{...}/
# - Per-run stdout/stderr -> run.log
# - Console: short START/DONE status only
# ============================================

try {
  chcp 65001 | Out-Null
  [Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
  $OutputEncoding = New-Object System.Text.UTF8Encoding($false)
} catch {}

$root = Split-Path -Parent $PSScriptRoot
$candidates = @(
  (Join-Path $root 'build\Release\morse.exe'),
  (Join-Path $root 'Release\morse.exe'),
  (Join-Path $root 'build\morse.exe')
)
$exe = $null
foreach ($c in $candidates) { if (Test-Path $c) { $exe = $c; break } }
if (-not $exe) { Write-Error "morse.exe not found"; exit 1 }
Write-Host "Using binary: $exe" -ForegroundColor Green

# PARAMETERS
$gammas = @(10, 20)
$Ks     = 4..8
$dts    = @('1e-4','1e-5','1e-6')
$N      = 2001     # grid points (adjust if needed)
$xmax   = 30       # domain max (adjust if needed)

$logMap = @{
  '1e-4' = 10000
  '1e-5' = 50000
  '1e-6' = 100000
}

# OUTPUT ROOTS
$resultsRoot = Join-Path $root 'results'
$baseTaylor  = Join-Path $resultsRoot 'taylor'
$baseMorse   = Join-Path $resultsRoot 'morse'
New-Item -ItemType Directory -Force -Path $baseTaylor | Out-Null
New-Item -ItemType Directory -Force -Path $baseMorse  | Out-Null

# PARALLELISM
$cores   = [Environment]::ProcessorCount
$maxJobs = [Math]::Max(1, [Math]::Min($cores - 1, 8))
Write-Host "Max parallel jobs: $maxJobs of $cores cores" -ForegroundColor Yellow

# Avoid oversubscription inside the binary/libs
$commonEnv = @{
  OMP_NUM_THREADS      = "1"
  MKL_NUM_THREADS      = "1"
  OPENBLAS_NUM_THREADS = "1"
}

# ------------------------------------------------------------------
# 1) PRECOMPUTE MORSE STATIC ONCE PER GAMMA (outside any task!)
# ------------------------------------------------------------------
foreach ($g in ($gammas | Select-Object -Unique)) {
  $gdir = Join-Path $baseMorse ("g{0}" -f $g)
  $potCsv = Join-Path $gdir 'morse_potential.csv'
  if (-not (Test-Path $potCsv)) {
    New-Item -ItemType Directory -Force -Path $gdir | Out-Null
    Write-Host ("[morse static] gamma={0} -> {1}" -f $g, $gdir) -ForegroundColor Yellow
    # run once to write morse_* files
    & $exe --gamma $g --N $N --xmax $xmax --outdir $gdir --stem ("morse_g{0}" -f $g) `
           *> (Join-Path $gdir 'run_static.log')
  } else {
    Write-Host ("[morse static] gamma={0} already exists, skipping" -f $g) -ForegroundColor DarkGray
  }
}

# ------------------------------------------------------------------
# 2) BUILD TASK LIST FOR TAYLOR EVOLUTION
# ------------------------------------------------------------------
$tasks = foreach ($g in $gammas) {
  foreach ($K in $Ks) {
    foreach ($dt in $dts) {
      [pscustomobject]@{ gamma=$g; K=$K; dt=$dt; log=$logMap[$dt] }
    }
  }
}

$active = @()
$evolveSwitch = '--evolve_only'  # matches your program's message "# --evolve_only: skipping ..."


function Start-OneTask {
  param($t, $exe, $baseTaylor, $envMap, $N, $xmax, $evolveSwitch)

  $gVal  = [int]$t.gamma
  $kVal  = [int]$t.K
  $dtVal = [string]$t.dt
  $logEv = [int]$t.log

  $dtTok = $dtVal.Replace('e-','e_').Replace('-','_')

  $sub    = Join-Path $baseTaylor ("g{0}\K{1}\dt_{2}" -f $gVal, $kVal, $dtTok)
  New-Item -ItemType Directory -Force -Path $sub | Out-Null

  $stem    = ("taylor_g{0}_K{1}_dt_{2}" -f $gVal, $kVal, $dtTok)
  $logFile = Join-Path $sub "run.log"

  Write-Host ("START [{0}] g={1} K={2} dt={3}" -f $stem, $gVal, $kVal, $dtVal) -ForegroundColor Cyan

  $sb = {
    param($exe,$gVal,$kVal,$dtVal,$logEv,$sub,$stem,$logFile,$envMap,$N,$xmax,$evolveSwitch)
    foreach ($k in $envMap.Keys) { [Environment]::SetEnvironmentVariable($k, $envMap[$k], "Process") }

    $cmd = @(
      $evolveSwitch, '--evolve','taylor',
      '--gamma', $gVal, '--K', $kVal, '--dt', $dtVal,
      '--tmax', 10, '--log-every', $logEv,
      '--N', $N, '--xmax', $xmax,
      '--outdir', $sub, '--stem', $stem
    ) -join ' '

    # Запишем команду в лог для проверки
    "[CMD] $exe $cmd" | Out-File -Encoding UTF8 -FilePath $logFile

    $start = Get-Date
    & $exe @($evolveSwitch,'--evolve','taylor','--gamma',$gVal,'--K',$kVal,'--dt',$dtVal,
             '--tmax',10,'--log-every',$logEv,'--N',$N,'--xmax',$xmax,'--outdir',$sub,'--stem',$stem) `
             *> $logFile
    $dur = (Get-Date) - $start
    "DONE  [$stem] g=$gVal K=$kVal dt=$dtVal (elapsed {0})" -f ($dur.ToString("hh\:mm\:ss"))
  }

  Start-Job -ScriptBlock $sb -ArgumentList $exe,$gVal,$kVal,$dtVal,$logEv,$sub,$stem,$logFile,$envMap,$N,$xmax,$evolveSwitch
}

# ...ниже в планировщике без изменений:
# $active += Start-OneTask -t $t -exe $exe -baseTaylor $baseTaylor -envMap $commonEnv -N $N -xmax $xmax -evolveSwitch $evolveSwitch


foreach ($t in $tasks) {
  while (($active | Where-Object State -eq 'Running').Count -ge $maxJobs) {
    Start-Sleep -Seconds 1
    $finished = $active | Where-Object State -ne 'Running'
    foreach ($j in $finished) {
      Receive-Job $j | ForEach-Object { Write-Host $_ -ForegroundColor Green }
      Remove-Job $j
      $active = $active | Where-Object Id -ne $j.Id
    }
  }
  $active += Start-OneTask -t $t -exe $exe -baseTaylor $baseTaylor -envMap $commonEnv -N $N -xmax $xmax -evolveSwitch $evolveSwitch
}

Wait-Job $active | Out-Null
foreach ($j in $active) {
  Receive-Job $j | ForEach-Object { Write-Host $_ -ForegroundColor Green }
  Remove-Job $j
}

Write-Host "All tasks finished. Output root: $resultsRoot" -ForegroundColor Green
