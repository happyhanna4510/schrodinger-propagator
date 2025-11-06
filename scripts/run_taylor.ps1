# ============================================
# Taylor runs in parallel (WinPS 5.1 & PS7 OK)
# - Static Morse written once per gamma -> results/morse/g{gamma}/
# - Taylor results -> results/taylor/g{gamma}/K{K}/dt_{...}/
# - Per-run stdout/stderr -> run.log (с [CMD] заголовком)
# - Для каждого прогона ~100 строк лога (auto --log-every)
# - tmax = 1 (сравнимо по физическому времени)
# ============================================
# ============================================
# Taylor runs in parallel (WinPS 5.1 & PS7 OK)
# - Static Morse once per gamma -> results/morse/g{gamma}/
# - Taylor results -> results/taylor/g{gamma}/K{K}/dt_{...}/
# - Per-run stdout/stderr -> run.log (with [CMD] header)
# - ~100 log lines per run (auto log_every)
# - tmax = 1
# ============================================

try {
  chcp 65001 | Out-Null
  [Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
  $OutputEncoding           = New-Object System.Text.UTF8Encoding($false)
} catch {}

$ErrorActionPreference = "Stop"

# --- locate exe (portable) ---
$repo = Resolve-Path (Join-Path $PSScriptRoot "..")
$exeCandidates = @(
  (Join-Path $repo "build\Release\morse.exe"),
  (Join-Path $repo "build\Debug\morse.exe"),
  (Join-Path $repo "build\morse.exe"),
  (Join-Path $repo "morse.exe")
) + (
  Get-ChildItem -Path $repo -Recurse -ErrorAction SilentlyContinue -Filter "morse*.exe" |
  Sort-Object LastWriteTime -Descending | Select-Object -First 1 | ForEach-Object { $_.FullName }
)
$exe = $exeCandidates | Where-Object { $_ -and (Test-Path $_) } | Select-Object -First 1
if (-not $exe) { throw "morse.exe not found. Build first." }
Write-Host "Using binary: $exe" -ForegroundColor Green

# --- PARAMETERS ---
$gammas = @(10, 20)
$Ks     = 4..8
$dts    = @('1e-6','1e-5','1e-4','1e-3','1e-2','1e-1')

$N      = 2001
$xmax   = 30
$tmax   = 1
$TARGET_LOG_LINES = 100

# --- OUTPUT ROOTS ---
$resultsRoot = Join-Path $repo 'results'
$baseTaylor  = Join-Path $resultsRoot 'taylor'
$baseMorse   = Join-Path $resultsRoot 'morse'
New-Item -ItemType Directory -Force -Path $baseTaylor | Out-Null
New-Item -ItemType Directory -Force -Path $baseMorse  | Out-Null

# --- PARALLELISM ---
$cores   = [Environment]::ProcessorCount
$maxJobs = [Math]::Max(1, [Math]::Min($cores - 1, 8))
Write-Host "Max parallel jobs: $maxJobs of $cores cores" -ForegroundColor Yellow

# --- avoid oversubscription ---
$commonEnv = @{
  OMP_NUM_THREADS      = "1"
  MKL_NUM_THREADS      = "1"
  OPENBLAS_NUM_THREADS = "1"
}

# --- 1) PRECOMPUTE MORSE STATIC ONCE PER GAMMA ---
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

# --- 2) BUILD TASK LIST (gamma × K × dt) ---
$tasks = foreach ($g in $gammas) {
  foreach ($K in $Ks) {
    foreach ($dt in $dts) {
      # nsteps and ~100 lines target
      $nsDbl  = [double]$tmax / [double]$dt
      $nsteps = [long]([math]::Max(1, [math]::Round($nsDbl)))
      $logEv  = [int]([math]::Max(1, [math]::Floor($nsteps / $TARGET_LOG_LINES)))
      [pscustomobject]@{ gamma=$g; K=$K; dt=$dt; log=$logEv }
    }
  }
}

$active = @()
$evolveSwitch = '--evolve_only'   # как и в RK4

function Start-TaylorTask {
  param($t, $exe, $baseTaylor, $envMap, $N, $xmax, $tmax, $evolveSwitch)

  $gVal  = [int]$t.gamma
  $kVal  = [int]$t.K
  $dtVal = [string]$t.dt
  $logEv = [int]$t.log

  $dtTok = ($dtVal -replace 'e-','e_') -replace '-','_'
  $sub   = Join-Path $baseTaylor ("g{0}\K{1}\dt_{2}" -f $gVal, $kVal, $dtTok)
  New-Item -ItemType Directory -Force -Path $sub | Out-Null

  $stem    = ("taylor_g{0}_K{1}_dt_{2}" -f $gVal, $kVal, $dtTok)
  $logFile = Join-Path $sub "run.log"

  Write-Host ("START [{0}] g={1} K={2} dt={3}" -f $stem, $gVal, $kVal, $dtVal) -ForegroundColor Cyan

  $sb = {
    param($exe,$gVal,$kVal,$dtVal,$logEv,$sub,$stem,$logFile,$envMap,$N,$xmax,$tmax,$evolveSwitch)

    foreach ($k in $envMap.Keys) { [Environment]::SetEnvironmentVariable($k, $envMap[$k], "Process") }

    $cmd = @(
      $evolveSwitch,'--evolve','taylor',
      '--gamma',$gVal,'--K',$kVal,'--dt',$dtVal,
      '--tmax',$tmax,'--log-every',$logEv,
      '--N',$N,'--xmax',$xmax,'--outdir',$sub,'--stem',$stem
    ) -join ' '

    # header
    Set-Content -Encoding UTF8 -Path $logFile -Value "[CMD] $exe $cmd"

    $start = Get-Date
    & $exe @($evolveSwitch,'--evolve','taylor',
             '--gamma',$gVal,'--K',$kVal,'--dt',$dtVal,
             '--tmax',$tmax,'--log-every',$logEv,
             '--N',$N,'--xmax',$xmax,'--outdir',$sub,'--stem',$stem) *>> $logFile
    $dur = (Get-Date) - $start

    # (опц.) Проваливай джобу, если бинарь вернул неноль
    if ($LASTEXITCODE -ne 0) {
      "EXIT CODE: $LASTEXITCODE" | Out-File -FilePath $logFile -Append -Encoding UTF8
      throw "Run failed: g=$gVal K=$kVal dt=$dtVal (exit $LASTEXITCODE)"
    }

    "DONE  [$stem] g=$gVal K=$kVal dt=$dtVal (elapsed {0})" -f ($dur.ToString("hh\:mm\:ss"))
  }

  Start-Job -ScriptBlock $sb -ArgumentList $exe,$gVal,$kVal,$dtVal,$logEv,$sub,$stem,$logFile,$envMap,$N,$xmax,$tmax,$evolveSwitch
}

# --- 3) SCHEDULER with throttling ---
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
  $active += Start-TaylorTask -t $t -exe $exe -baseTaylor $baseTaylor -envMap $commonEnv -N $N -xmax $xmax -tmax $tmax -evolveSwitch $evolveSwitch
}

Wait-Job $active | Out-Null
foreach ($j in $active) {
  Receive-Job $j | ForEach-Object { Write-Host $_ -ForegroundColor Green }
  Remove-Job $j
}

Write-Host "All Taylor tasks finished. Output root: $resultsRoot" -ForegroundColor Green
