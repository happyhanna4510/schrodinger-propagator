# =============================
# run_field.ps1  (batch launcher for morse.exe)
# =============================

$exe         = ".\build\Release\morse.exe"
$resultsRoot = ".\results"

# --- fixed numerical params ---
$N    = 2001
$xmax = 20
$tmax = 10
$k0   = 0
$logEvery = 100  

# --- methods ---
$methods = @("cheb")
# $methods = @("cheb", "rk4", "taylor")

$taylorK = 4


$U0List    = @(-1, 1, 5, -5, -100, -30)
$gammaList = @(10)
$dtList    = @("1e-3")
# $dtList    = @("1e-5","1e-4","1e-3")

# =============================

function Ensure-Dir($path) {
    if (-not (Test-Path -LiteralPath $path)) {
        New-Item -ItemType Directory -Path $path | Out-Null
    }
}

Write-Host "=== Start batch run ==="

foreach ($U0 in $U0List) {

    $u0FolderName = "U0_$U0"
    $u0FolderPath = Join-Path $resultsRoot $u0FolderName
    Ensure-Dir $u0FolderPath

    foreach ($gamma in $gammaList) {

        $gFolderName = "g$gamma"
        $gFolderPath = Join-Path $u0FolderPath $gFolderName
        Ensure-Dir $gFolderPath

        foreach ($dtStr in $dtList) {

            $dtVal        = [double]$dtStr
            $dtFolderName = "dt$dtStr"
            $dtFolderPath = Join-Path $gFolderPath $dtFolderName
            Ensure-Dir $dtFolderPath

            Write-Host ""
            Write-Host ">>> U0 = $U0, gamma = $gamma, dt = $dtStr"
            Write-Host "    outdir = $dtFolderPath"
            Write-Host "    log-every = $logEvery steps"

            foreach ($method in $methods) {

                # method tag for filename
                if ($method -eq "taylor") {
                    $methodTag = "taylor_K$($taylorK)"
                } else {
                    $methodTag = $method
                }

                # CSV filename (same convention as before)
                $csvName = "{0}_dt{1}_g{2}_N{3}_x{4}_U0_{5}_k0_{6}_tmax_{7}.csv" -f `
                    $methodTag, $dtStr, $gamma, $N, $xmax, $U0, $k0, $tmax

                $csvPath = Join-Path $dtFolderPath $csvName

                # common args (without --evolve / --K)
                $commonArgs = @(
                    "--gamma",     $gamma,
                    "--N",         $N,
                    "--xmax",      $xmax,
                    "--tmax",      $tmax,
                    "--dt",        $dtVal,
                    "--init",      "complex-gauss",
                    "--wide",
                    "--k0",        $k0,
                    "--U0",        $U0,
                    "--outdir",    $dtFolderPath,
                    "--csv",       $csvName,
                    "--log-every", $logEvery    # <=== ВОТ ЗДЕСЬ МЫ ПРОТАСКИВАЕМ ФЛАГ
                )

                if ($method -eq "taylor") {
                    $args = @(
                        "--evolve", "taylor",
                        "--K",      $taylorK
                    ) + $commonArgs
                } else {
                    $args = @(
                        "--evolve", $method
                    ) + $commonArgs
                }

                Write-Host "  -> $methodTag  (CSV: $csvName)"
                & $exe @args
                if ($LASTEXITCODE -ne 0) {
                    Write-Warning ("    !!! morse.exe exited with code {0}" -f $LASTEXITCODE)
                }
            }
        }
    }
}

Write-Host "=== Done ==="
