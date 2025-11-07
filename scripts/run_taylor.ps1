# ============================================
# Taylor runs in parallel (WinPS 5.1 & PS7 OK)
# - Static Morse once per gamma -> results/morse/g{gamma}/
# - Taylor -> results/taylor/g{gamma}/K{K}/dt_{...}/
# - Per-run stdout/stderr -> run.log (with [CMD] header)
# - CSV -> unique long filename per dt
# - ~100 log lines per run (auto --log-every), tmax = 1
# - Parallel status banners + per-job elapsed time
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

# Complex Gaussian defaults (for unique file tags)
$Init  = 'complex-gauss'
$X0    = 0
$Sigma = 1
$K0    = 10
$U0    = 0
$PktTag = ("cg_x0_{0}_s{1}_k{2}_u{3}" -f $X0, $Sigma, $K0, $U0).Replace('.', 'p')

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

# --- avoid oversubscription in libs ---
$commonEnv = @{
  OMP_NUM_THREADS      = "1"
  MKL_NUM_THREADS      = "1"
  OPENBLAS_NUM_THREADS = "1"
}

# --- 1) PRECOMPUTE MORSE STATIC ONCE PER GAMMA (NO --stem) ---
foreach ($g in ($gammas | Select-Object -Unique)) {
  $gdir   = Join-Path $baseMorse ("g{0}" -f $g)
  $potCsv = Join-Path $gdir 'morse_potential.csv'
  if (-not (Test-Path $gdir)) {
    New-Item -ItemType Directory -Force -Path $gdir | Out-Null
  }

  if (-not (Test-Path $potCsv)) {
    Write-Host ("[morse static] gamma={0} -> {1}" -f $g, $gdir) -ForegroundColor Yellow
    & $exe --gamma $g --N $N --xmax $xmax --outdir $gdir `
           *> (Join-Path $gdir 'run_static.log')
  } else {
    Write-Host ("[morse static] gamma={0} already exists, skipping" -f $g) -ForegroundColor DarkGray
  }
}

# --- 2) BUILD TASK LIST (gamma × K × dt) ---
$tasks = foreach ($g in $gammas) {
  foreach ($K in $Ks) {
    foreach ($dt in $dts) {
      $nsDbl  = [double]$tmax / [double]$dt
      $nsteps = [long]([math]::Max(1, [math]::Round($nsDbl)))
      $logEv  = [int]([math]::Max(1, [math]::Floor($nsteps / $TARGET_LOG_LINES)))
      [pscustomobject]@{ gamma=$g; K=$K; dt=$dt; log=$logEv }
    }
  }
}
$total = $tasks.Count
Write-Host ("[taylor] launching {0} tasks with up to {1} parallel jobs..." -f $total, $maxJobs) -ForegroundColor Yellow

# --- Start-TaylorTask returns a Job context ---
function Start-TaylorTask {
  param($t, $exe, $baseTaylor, $envMap, $N, $xmax, $tmax, $evolveSwitch,
        $Init,$X0,$Sigma,$K0,$U0,$PktTag)

  $gVal  = [int]$t.gamma
  $kVal  = [int]$t.K
  $dtVal = [string]$t.dt
  $logEv = [int]$t.log

  $dtTok = ($dtVal -replace 'e-','e_') -replace '-','_'

  # results/taylor/g{g}/K{k}/dt_{dtTok}/
  $sub = Join-Path $baseTaylor ("g{0}\K{1}\dt_{2}" -f $gVal, $kVal, $dtTok)
  New-Item -ItemType Directory -Force -Path $sub | Out-Null

  # unique long CSV name
  $csvBase = ("taylor_g{0}_K{1}_dt_{2}_N{3}_xmax{4}_tmax{5}_{6}" -f `
              $gVal,$kVal,$dtTok,$N,($xmax.ToString().Replace('.','p')),$tmax,$PktTag)

  # temp outdir inside dt_
  $tmpOut = Join-Path $sub ("__tmp_" + [guid]::NewGuid().ToString("N"))
  New-Item -ItemType Directory -Force -Path $tmpOut | Out-Null

  # per-run log in dt_ folder
  $logFile = Join-Path $sub "run.log"

  Write-Host ("START [taylor] g={0}  K={1}  dt={2}  -> {3}.csv" -f $gVal,$kVal,$dtVal,$csvBase) -ForegroundColor Cyan

  $sb = {
    param($exe,$gVal,$kVal,$dtVal,$logEv,$tmpOut,$logFile,$envMap,$N,$xmax,$tmax,$evolveSwitch,$Init,$X0,$Sigma,$K0,$U0)

    foreach ($k in $envMap.Keys) { [Environment]::SetEnvironmentVariable($k, $envMap[$k], "Process") }

    $cmd = @(
      $evolveSwitch,'--evolve','taylor',
      '--gamma',$gVal,'--K',$kVal,'--dt',$dtVal,
      '--tmax',$tmax,'--log-every',$logEv,
      '--N',$N,'--xmax',$xmax,
      '--outdir',$tmpOut,
      '--init',$Init,'--x0',$X0,'--sigma',$Sigma,'--k0',$K0,'--U0',$U0
    ) -join ' '

    Set-Content -Encoding UTF8 -Path $logFile -Value ("[CMD] {0} {1}" -f $exe, $cmd)

    $start = Get-Date
    & $exe @(
      $evolveSwitch,'--evolve','taylor',
      '--gamma',$gVal,'--K',$kVal,'--dt',$dtVal,
      '--tmax',$tmax,'--log-every',$logEv,
      '--N',$N,'--xmax',$xmax,
      '--outdir',$tmpOut,
      '--init',$Init,'--x0',$X0,'--sigma',$Sigma,'--k0',$K0,'--U0',$U0
    ) *>> $logFile
    $dur = (Get-Date) - $start

    if ($LASTEXITCODE -ne 0) {
      "EXIT CODE: $LASTEXITCODE" | Out-File -FilePath $logFile -Append -Encoding UTF8
      throw ("Run failed: g={0} K={1} dt={2} (exit {3})" -f $gVal,$kVal,$dtVal,$LASTEXITCODE)
    }

    # return metadata (captured by Receive-Job)
    [pscustomobject]@{ Tmp=$tmpOut; Took=$dur.ToString("hh\:mm\:ss") }
  }

  $job = Start-Job -ScriptBlock $sb -ArgumentList `
    $exe,$gVal,$kVal,$dtVal,$logEv,$tmpOut,$logFile,$commonEnv,$N,$xmax,$tmax,'--evolve_only',$Init,$X0,$Sigma,$K0,$U0

  return [pscustomobject]@{
    Job     = $job
    TmpOut  = $tmpOut
    SubDir  = $sub
    CsvBase = $csvBase
    G       = $gVal
    K       = $kVal
    DT      = $dtVal
  }
}

function Finalize-Job {
  param($ctx, $took)

  # pick the largest CSV as the primary result
  $primary = Get-ChildItem -Path $ctx.TmpOut -Filter "*.csv" -ErrorAction SilentlyContinue |
             Sort-Object Length -Descending | Select-Object -First 1
  if ($primary) {
    $targetCsv = Join-Path $ctx.SubDir ("{0}.csv" -f $ctx.CsvBase)
    Move-Item -Force $primary.FullName $targetCsv
  }
  # optional: move other CSVs with suffix
  # Get-ChildItem -Path $ctx.TmpOut -Filter "*.csv" |
  #   Where-Object { $_.FullName -ne $primary.FullName } |
  #   ForEach-Object {
  #     $tgt = Join-Path $ctx.SubDir ("{0}__{1}" -f $ctx.CsvBase, $_.Name)
  #     Move-Item -Force $_.FullName $tgt
  #   }

  Remove-Item -Recurse -Force $ctx.TmpOut
  if (-not $took) { $took = "unknown" }
  Write-Host ("DONE  [taylor] g={0}  K={1}  dt={2}  -> {3}.csv  (elapsed {4})" -f `
              $ctx.G,$ctx.K,$ctx.DT,$ctx.CsvBase,$took) -ForegroundColor Green
}

# --- 3) SCHEDULER with throttling (parallel) ---
$active   = @()
$started  = 0
$finished = 0

function Show-ParallelStatus {
  param($started,$finished,$total,$running)
  Write-Host ("[parallel] started={0}/{1}; finished={2}/{1}; running={3}" -f `
              $started,$total,$finished,$running) -ForegroundColor DarkCyan
}

foreach ($t in $tasks) {
  while (($active | Where-Object { $_.Job.State -eq 'Running' }).Count -ge $maxJobs) {
    Start-Sleep -Seconds 1
    $done = $active | Where-Object { $_.Job.State -ne 'Running' }
    foreach ($ctx in $done) {
      $out  = Receive-Job $ctx.Job
      $took = ($out | Where-Object { $_.PSObject.Properties.Name -contains 'Took' } | Select-Object -First 1).Took
      Finalize-Job $ctx $took
      Remove-Job $ctx.Job
      $finished++
      $active = $active | Where-Object { $_.Job.Id -ne $ctx.Job.Id }
      Show-ParallelStatus -started $started -finished $finished -total $total -running ($active.Count)
    }
  }

  $ctx = Start-TaylorTask -t $t `
    -exe $exe -baseTaylor $baseTaylor -envMap $commonEnv `
    -N $N -xmax $xmax -tmax $tmax -evolveSwitch '--evolve_only' `
    -Init $Init -X0 $X0 -Sigma $Sigma -K0 $K0 -U0 $U0 -PktTag $PktTag

  $active += $ctx
  $started++
  Show-ParallelStatus -started $started -finished $finished -total $total -running ($active.Count)
}

# drain remaining jobs
Wait-Job ($active | ForEach-Object { $_.Job }) | Out-Null
foreach ($ctx in @($active)) {
  $out  = Receive-Job $ctx.Job
  $took = ($out | Where-Object { $_.PSObject.Properties.Name -contains 'Took' } | Select-Object -First 1).Took
  Finalize-Job $ctx $took
  Remove-Job $ctx.Job
  $finished++
  $active = $active | Where-Object { $_.Job.Id -ne $ctx.Job.Id }
  Show-ParallelStatus -started $started -finished $finished -total $total -running ($active.Count)
}

Write-Host ("All Taylor tasks finished. Total: {0}. Max parallel: {1}. Output: {2}" -f $total, $maxJobs, $resultsRoot) -ForegroundColor Green
