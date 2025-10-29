# ============================================
# Chebyshev runs in parallel (WinPS 5.1 & PS7 OK)
# - Static Morse once per gamma -> results/morse/g{gamma}/
# - Cheb results -> results/cheb/g{gamma}/dt_{...}/tol_{...}/
# - Per-run stdout/stderr -> run.log (with [CMD] header)
# ============================================

try {
  chcp 65001 | Out-Null
  [Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
  $OutputEncoding = New-Object System.Text.UTF8Encoding($false)
} catch {}

# locate exe
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
$dts    = @('1e-4','1e-5','1e-6')
$tols   = @('1e-11')
$N      = 2001   # adjust if needed
$xmax   = 30     # adjust if needed
$tmax   = 10

# log_every by dt
$logMap = @{
  '1e-4' = 10000
  '1e-5' = 50000
  '1e-6' = 100000
}

# OUTPUT ROOTS
$resultsRoot = Join-Path $root 'results'
$baseCheb    = Join-Path $resultsRoot 'cheb'
$baseMorse   = Join-Path $resultsRoot 'morse'
New-Item -ItemType Directory -Force -Path $baseCheb  | Out-Null
New-Item -ItemType Directory -Force -Path $baseMorse | Out-Null

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

# 1) PRECOMPUTE MORSE STATIC ONCE PER GAMMA
foreach ($g in ($gammas | Select-Object -Unique)) {
  $gdir   = Join-Path $baseMorse ("g{0}" -f $g)
  $potCsv = Join-Path $gdir 'morse_potential.csv'
  if (-not (Test-Path $potCsv)) {
    New-Item -ItemType Directory -Force -Path $gdir | Out-Null
    Write-Host ("[morse static] gamma={0} -> {1}" -f $g, $gdir) -ForegroundColor Yellow
    & $exe --gamma $g --N $N --xmax $xmax --outdir $gdir --stem ("morse_g{0}" -f $g) `
           *> (Join-Path $gdir 'run_static.log')
  } else {
    Write-Host ("[morse static] gamma={0} already exists, skipping" -f $g) -ForegroundColor DarkGray
  }
}

# 2) BUILD TASKS
$tasks = foreach ($g in $gammas) {
  foreach ($dt in $dts) {
    foreach ($tol in $tols) {
      [pscustomobject]@{ gamma=$g; dt=$dt; tol=$tol; log=$logMap[$dt] }
    }
  }
}

$active = @()
$evolveSwitch = '--evolve_only'  # your binary prints "# --evolve_only: skipping ..."

function Start-ChebTask {
  param($t, $exe, $baseCheb, $envMap, $N, $xmax, $tmax, $evolveSwitch)

  $gVal   = [int]$t.gamma
  $dtVal  = [string]$t.dt
  $tolVal = [string]$t.tol
  $logEv  = [int]$t.log

  $dtTok  = $dtVal.Replace('e-','e_').Replace('-','_')
  $tolTok = $tolVal.Replace('e-','e_').Replace('-','_')  # e.g., 1e-11 -> 1e_11

  $sub = Join-Path $baseCheb ("g{0}\dt_{1}\tol_{2}" -f $gVal, $dtTok, $tolTok)
  New-Item -ItemType Directory -Force -Path $sub | Out-Null

  $stem    = ("cheb_g{0}_dt_{1}_tol_{2}" -f $gVal, $dtTok, $tolTok)
  $logFile = Join-Path $sub "run.log"

  Write-Host ("START [{0}] g={1} dt={2} tol={3}" -f $stem, $gVal, $dtVal, $tolVal) -ForegroundColor Cyan

  $sb = {
    param($exe,$gVal,$dtVal,$tolVal,$logEv,$sub,$stem,$logFile,$envMap,$N,$xmax,$tmax,$evolveSwitch)

    foreach ($k in $envMap.Keys) { [Environment]::SetEnvironmentVariable($k, $envMap[$k], "Process") }

    $cmd = @(
      $evolveSwitch,'--evolve','cheb',
      '--gamma',$gVal,'--dt',$dtVal,'--tmax',$tmax,
      '--tol',$tolVal,'--K',0,'--log-every',$logEv,
      '--N',$N,'--xmax',$xmax,'--outdir',$sub,'--stem',$stem
    ) -join ' '

    # write command header; then append process output
    Set-Content -Encoding UTF8 -Path $logFile -Value "[CMD] $exe $cmd"

    $start = Get-Date
    & $exe @($evolveSwitch,'--evolve','cheb','--gamma',$gVal,'--dt',$dtVal,'--tmax',$tmax,
             '--tol',$tolVal,'--K',0,'--log-every',$logEv,'--N',$N,'--xmax',$xmax,
             '--outdir',$sub,'--stem',$stem) *>> $logFile
    $dur = (Get-Date) - $start

    "DONE  [$stem] g=$gVal dt=$dtVal tol=$tolVal (elapsed {0})" -f ($dur.ToString("hh\:mm\:ss"))
  }

  Start-Job -ScriptBlock $sb -ArgumentList $exe,$gVal,$dtVal,$tolVal,$logEv,$sub,$stem,$logFile,$envMap,$N,$xmax,$tmax,$evolveSwitch
}

# scheduler with throttling
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
  $active += Start-ChebTask -t $t -exe $exe -baseCheb $baseCheb -envMap $commonEnv -N $N -xmax $xmax -tmax $tmax -evolveSwitch $evolveSwitch
}

Wait-Job $active | Out-Null
foreach ($j in $active) {
  Receive-Job $j | ForEach-Object { Write-Host $_ -ForegroundColor Green }
  Remove-Job $j
}

Write-Host "All Chebyshev tasks finished. Output root: $resultsRoot" -ForegroundColor Green
