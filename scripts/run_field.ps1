# =============================
# run_field_parallel.ps1
# =============================

$exe         = ".\build\Release\morse.exe"
$resultsRoot = ".\results"

# --- fixed numerical params ---
$N        = 2001
$xmax     = 20
$tmax     = 100
$k0       = 0
$logEvery = 1000  

# --- methods ---
$methods  = @("cheb", "rk4", "taylor")
$taylorK  = 4

$U0List    = @(-1, 1)
$gammaList = @(10)
$dtList    = @("1e-5","1e-4","1e-3")

$throttleLimit = 7   # использовать 7 потоков из 8

function Ensure-Dir($path) {
    [System.IO.Directory]::CreateDirectory($path) | Out-Null
}

Write-Host "=== Build job list ==="

$jobs = @()

foreach ($U0 in $U0List) {
    foreach ($gamma in $gammaList) {
        foreach ($dtStr in $dtList) {
            foreach ($method in $methods) {

                $jobs += [PSCustomObject]@{
                    U0     = $U0
                    Gamma  = $gamma
                    DtStr  = $dtStr
                    Method = $method
                }
            }
        }
    }
}

Write-Host "Total jobs: $($jobs.Count)"

$jobs | ForEach-Object -Parallel {

    $U0     = $_.U0
    $gamma  = $_.Gamma
    $dtStr  = $_.DtStr
    $method = $_.Method

    $N        = $using:N
    $xmax     = $using:xmax
    $tmax     = $using:tmax
    $k0       = $using:k0
    $logEvery = $using:logEvery
    $taylorK  = $using:taylorK

    $resultsRoot = $using:resultsRoot
    $exe         = $using:exe

    function Ensure-Dir($path) {
        [System.IO.Directory]::CreateDirectory($path) | Out-Null
    }

    $u0FolderName = "U0_$U0"
    $u0FolderPath = Join-Path $resultsRoot $u0FolderName
    Ensure-Dir $u0FolderPath

    $gFolderName = "g$gamma"
    $gFolderPath = Join-Path $u0FolderPath $gFolderName
    Ensure-Dir $gFolderPath

    $dtFolderName = "dt$dtStr"
    $dtFolderPath = Join-Path $gFolderPath $dtFolderName
    Ensure-Dir $dtFolderPath

    $dtVal = [double]$dtStr

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
        "--outdir",    $dtFolderPath,
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

    Write-Host "[$([System.Threading.Thread]::CurrentThread.ManagedThreadId)] U0=$U0 gamma=$gamma dt=$dtStr method=$methodTag"
    & $exe @args
    if ($LASTEXITCODE -ne 0) {
        Write-Warning ("    !!! morse.exe exited with code {0} (U0={1}, g={2}, dt={3}, method={4})" -f `
            $LASTEXITCODE, $U0, $gamma, $dtStr, $methodTag)
    }

} -ThrottleLimit $throttleLimit

Write-Host "=== All jobs finished ==="
