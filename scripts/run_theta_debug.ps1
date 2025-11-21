# =============================
# run_theta_debug.ps1
# =============================
# Запускает morse.exe с включённой диагностикой фазы (THETA_DEBUG)
# для выбранных методов и значений U0.
# Сохраняет полный консольный вывод в отдельные .log файлы,
# чтобы потом удобно парсить их в Python.

$exe         = ".\build\Release\morse.exe"
$resultsRoot = ".\results"

# --- численные параметры (как в обычном запуске) ---
$N    = 2001
$xmax = 20
$tmax = 10
$dtStr = "1e-5"
$dtVal = [double]$dtStr
$gamma = 10
$k0    = 0

# --- методы, которые хотим проверить ---
$methods = @("cheb", "rk4")   # можно убрать/добавить по желанию
#$methods = @("rk4", "taylor")   # можно убрать/добавить по желанию

$taylorK = 4

# --- список U0 ---
# Сейчас можно оставить только -1, позже расширить:
# $U0List = @(-1, 1, 5, -5, -30, -100)
$U0List = @(-100,-1,1)

# --- настройки отладки фазы ---
$thetaDebugEnabled = "1"        # THETA_DEBUG=1 включает диагностику
$thetaDebugWindow  = "3.5,4.5"  # окно по времени для подробных логов overlap

# Устанавливаем переменные окружения для текущей сессии PowerShell
$env:THETA_DEBUG        = $thetaDebugEnabled
$env:THETA_DEBUG_WINDOW = $thetaDebugWindow

function Ensure-Dir($path) {
    if (-not (Test-Path -LiteralPath $path)) {
        New-Item -ItemType Directory -Path $path | Out-Null
    }
}

Write-Host "=== Theta debug batch run ==="
Write-Host "THETA_DEBUG=$env:THETA_DEBUG"
Write-Host "THETA_DEBUG_WINDOW=$env:THETA_DEBUG_WINDOW"
Write-Host ""

foreach ($U0 in $U0List) {

    $u0FolderName = "U0_$U0"
    $u0FolderPath = Join-Path $resultsRoot $u0FolderName
    Ensure-Dir $u0FolderPath

    $gFolderName = "g$gamma"
    $gFolderPath = Join-Path $u0FolderPath $gFolderName
    Ensure-Dir $gFolderPath

    $dtFolderName = "dt$dtStr"
    $dtFolderPath = Join-Path $gFolderPath $dtFolderName
    Ensure-Dir $dtFolderPath

    Write-Host ">>> U0 = $U0, gamma = $gamma, dt = $dtStr"
    Write-Host "    outdir = $dtFolderPath"
    Write-Host ""

    foreach ($method in $methods) {

        if ($method -eq "taylor") {
            $methodTag = "taylor_K$($taylorK)"
        } else {
            $methodTag = $method
        }

        # Имя CSV (как у тебя в обычных запусках)
        $csvName = "{0}_dt{1}_g{2}_N{3}_x{4}_U0_{5}_k0_{6}_tmax_{7}.csv" -f `
            $methodTag, $dtStr, $gamma, $N, $xmax, $U0, $k0, $tmax

        $csvPath = Join-Path $dtFolderPath $csvName

        # Имя лог-файла с диагностикой theta
        $logName = "theta_debug_{0}_dt{1}_g{2}_U0_{3}.log" -f `
            $methodTag, $dtStr, $gamma, $U0
        $logPath = Join-Path $dtFolderPath $logName

        # Общие аргументы
        $commonArgs = @(
            "--gamma",  $gamma,
            "--N",      $N,
            "--xmax",   $xmax,
            "--tmax",   $tmax,
            "--dt",     $dtVal,
            "--init",   "complex-gauss",
            "--k0",     $k0,
            "--U0",     $U0,
            "--outdir", $dtFolderPath,
            "--csv",    $csvName
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

        Write-Host "  -> method = $methodTag"
        Write-Host "     CSV : $csvName"
        Write-Host "     LOG : $logName"
        Write-Host ""

        # Запускаем morse.exe и ПЕРЕНАПРАВЛЯЕМ весь вывод в лог
        # *> перенаправляет все потоки (stdout, stderr и т.п.)
        & $exe @args *> $logPath

        if ($LASTEXITCODE -ne 0) {
            Write-Warning ("    !!! morse.exe exited with code {0}" -f $LASTEXITCODE)
        } else {
            Write-Host "    OK"
        }

        Write-Host ""
    }
}

Write-Host "=== Theta debug batch done ==="
