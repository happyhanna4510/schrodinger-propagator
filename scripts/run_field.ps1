# ============================================
# Parallel field runs for cheb / rk4 / taylor
# Layout:
#   results/
#     U0_<U0>/g<gamma>/dt<dt>/
#       run.log
#       <methodTag>_dt<dt>_g<gamma>_N<N>_x<xmax>_U0_<U0>_k0_<k0>_tmax_<tmax>.csv
#       (plus wide-files only for cheb)
# ============================================

try {
  chcp 65001 | Out-Null
  [Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
  $OutputEncoding           = New-Object System.Text.UTF8Encoding($false)
} catch {}

$ErrorActionPreference = "Stop"

# --- locate exe ---
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

# --- NUMERICAL PARAMS ---
$N      = 2001
$xmax   = 20
$tmax   = 100
$k0     = 0

# methods
$methods = @("cheb", "rk4", "taylor")
$taylorK = 4

# physics
$U0List    = @(-1, 1)
$gammaList = @(10)
$dtList    = @("1e-5","1e-4","1e-3")

$TARGET_LOG_LINES = 200

# --- OUTPUT ROOT ---
$resultsRoot = Join-Path $repo 'results'
New-Item -ItemType Directory -Force -Path $resultsRoot | Out-Null

# --- PARALLELISM ---
$cores   = [Environment]::ProcessorCount
$maxJobs = [Math]::Max(1, [Math]::Min($cores - 1, 7))   # <= ТУТ ТЕПЕРЬ 7
Write-Host "Max parallel jobs: $maxJobs of $cores cores" -ForegroundColor Yellow

# --- avoid oversubscription in BLAS ---
$commonEnv = @{
  OMP_NUM_THREADS      = "1"
  MKL_NUM_THREADS      = "1"
  OPENBLAS_NUM_THREADS = "1"
}

# --- build task list ---
$tasks = foreach ($U0 in $U0List) {
  foreach ($gamma in $gammaList) {
    foreach ($dtStr in $dtList) {
      $nsteps  = [double]$tmax / [double]$dtStr
      $logEv   = [int]([math]::Max(1, [math]::Floor($nsteps / $TARGET_LOG_LINES)))
      [pscustomobject]@{
        U0       = $U0
        Gamma    = $gamma
        DtStr    = $dtStr
        LogEvery = $logEv
      }
    }
  }
}

$total = $tasks.Count
Write-Host ("[field] launching {0} tasks (U0,gamma,dt combos) with up to {1} parallel jobs..." -f $total, $maxJobs) -ForegroundColor Yellow

function Start-FieldTask {
  param(
    $t, $exe, $resultsRoot, $methods,
    $N, $xmax, $tmax, $k0, $envMap, $taylorK
  )

  $U0       = [double]$t.U0
  $gamma    = [int]$t.Gamma
  $dtStr    = [string]$t.DtStr
  $logEvery = [int]$t.LogEvery
  $dtVal    = [double]$dtStr

  $u0Folder = Join-Path $resultsRoot ("U0_{0}" -f $U0)
  $gFolder  = Join-Path $u0Folder   ("g{0}"    -f $gamma)
  $dtFolder = Join-Path $gFolder    ("dt{0}"   -f $dtStr)

  New-Item -ItemType Directory -Force -Path $u0Folder | Out-Null
  New-Item -ItemType Directory -Force -Path $gFolder  | Out-Null
  New-Item -ItemType Directory -Force -Path $dtFolder | Out-Null

  $job = Start-Job -ScriptBlock {
    param(
      $exe, $U0, $gamma, $dtStr, $dtVal, $dtFolder,
      $methods, $N, $xmax, $tmax, $k0, $logEvery,
      $envMap, $taylorK
    )

    foreach ($k in $envMap.Keys) {
      [Environment]::SetEnvironmentVariable($k, $envMap[$k], "Process")
    }

    $runLog = Join-Path $dtFolder "run.log"
    Set-Content -Encoding UTF8 -Path $runLog -Value ("[FIELD] U0={0} gamma={1} dt={2} tmax={3}" -f $U0,$gamma,$dtStr,$tmax)

    $startTotal = Get-Date

    foreach ($method in $methods) {

      if ($method -eq "taylor") {
        $methodTag = "taylor_K$($taylorK)"
      } else {
        $methodTag = $method
      }

      $csvName = "{0}_dt{1}_g{2}_N{3}_x{4}_U0_{5}_k0_{6}_tmax_{7}.csv" -f `
                 $methodTag, $dtStr, $gamma, $N, $xmax, $U0, $k0, $tmax

      $baseArgs = @(
        "--gamma",     $gamma,
        "--N",         $N,
        "--xmax",      $xmax,
        "--tmax",      $tmax,
        "--dt",        $dtVal,
        "--init",      "complex-gauss",
        "--k0",        $k0,
        "--U0",        $U0,
        "--outdir",    $dtFolder,
        "--csv",       $csvName,
        "--log-every", $logEvery
      )

      switch ($method) {
        "cheb" {
          $args = @(
            "--evolve", "cheb",
            "--wide"
          ) + $baseArgs
        }
        "rk4" {
          $args = @(
            "--evolve", "rk4"
          ) + $baseArgs
        }
        "taylor" {
          $args = @(
            "--evolve", "taylor",
            "--K",      $taylorK
          ) + $baseArgs
        }
      }

      Add-Content -Encoding UTF8 -Path $runLog -Value ("[CMD] {0} {1}" -f $exe, ($args -join ' '))

      $start = Get-Date
      & $exe $args *>> $runLog
      $dur = (Get-Date) - $start

      if ($LASTEXITCODE -ne 0) {
        Add-Content -Encoding UTF8 -Path $runLog -Value ("[ERROR] exit code {0} for method {1}" -f $LASTEXITCODE,$methodTag)
        throw ("Field run failed: U0={0} gamma={1} dt={2} method={3} (exit {4})" -f $U0,$gamma,$dtStr,$methodTag,$LASTEXITCODE)
      } else {
        Add-Content -Encoding UTF8 -Path $runLog -Value ("[OK] {0} finished in {1}" -f $methodTag,$dur.ToString("hh\:mm\:ss"))
      }
    }

    $totalDur = (Get-Date) - $startTotal
    [pscustomobject]@{
      U0    = $U0
      Gamma = $gamma
      DtStr = $dtStr
      Took  = $totalDur.ToString("hh\:mm\:ss")
      Path  = $dtFolder
    }

  } -ArgumentList $exe,$U0,$gamma,$dtStr,$dtVal,$dtFolder,`
                 $methods,$N,$xmax,$tmax,$k0,$logEvery,`
                 $commonEnv,$taylorK

  return [pscustomobject]@{
    Job     = $job
    U0      = $U0
    Gamma   = $gamma
    DtStr   = $dtStr
    DtPath  = $dtFolder
  }
}

function Show-ParallelStatus {
  param($started,$finished,$total,$running)
  Write-Host ("[parallel] started={0}/{1}; finished={2}/{1}; running={3}" -f `
              $started,$total,$finished,$running) -ForegroundColor DarkCyan
}

$active   = @()
$started  = 0
$finished = 0

foreach ($t in $tasks) {

  while (($active | Where-Object { $_.Job.State -eq 'Running' }).Count -ge $maxJobs) {
    Start-Sleep -Seconds 1
    $done = $active | Where-Object { $_.Job.State -ne 'Running' }
    foreach ($ctx in $done) {
      $out = Receive-Job $ctx.Job
      if ($out) {
        $took = $out.Took
        Write-Host ("DONE  [field] U0={0} gamma={1} dt={2} -> {3} (elapsed {4})" -f `
                    $out.U0,$out.Gamma,$out.DtStr,$out.Path,$took) -ForegroundColor Green
      }
      Remove-Job $ctx.Job
      $finished++
      $active = $active | Where-Object { $_.Job.Id -ne $ctx.Job.Id }
      Show-ParallelStatus -started $started -finished $finished -total $total -running ($active.Count)
    }
  }

  $ctx = Start-FieldTask -t $t `
    -exe $exe -resultsRoot $resultsRoot -methods $methods `
    -N $N -xmax $xmax -tmax $tmax -k0 $k0 -envMap $commonEnv -taylorK $taylorK

  $active += $ctx
  $started++
  Show-ParallelStatus -started $started -finished $finished -total $total -running ($active.Count)
}

Wait-Job ($active | ForEach-Object { $_.Job }) | Out-Null
foreach ($ctx in @($active)) {
  $out = Receive-Job $ctx.Job
  if ($out) {
    $took = $out.Took
    Write-Host ("DONE  [field] U0={0} gamma={1} dt={2} -> {3} (elapsed {4})" -f `
                $out.U0,$out.Gamma,$out.DtStr,$out.Path,$took) -ForegroundColor Green
  }
  Remove-Job $ctx.Job
  $finished++
  $active = $active | Where-Object { $_.Job.Id -ne $ctx.Job.Id }
  Show-ParallelStatus -started $started -finished $finished -total $total -running ($active.Count)
}

Write-Host ("All field tasks finished. Total: {0}. Max parallel: {1}. Output root: {2}" -f $total, $maxJobs, $resultsRoot) -ForegroundColor Green


