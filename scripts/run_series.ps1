# Путь к exe
$exe = ".\build\Release\morse.exe"

# Фиксированные параметры
$gamma    = 10
$N        = 2001
$xmax     = 30
$Umax     = 0.1
$tmax     = 10.0
$logEvery = 10000    # реже = быстрее

# Сетки значений
$dts = @("1e-6","1e-5","1e-4")
$Ks  = 4..8

#New-Item -ItemType Directory -Force -Path "results" | Out-Null

foreach ($K in $Ks) {
    foreach ($dt in $dts) {
        if ($dt -eq "1e-4") { $logEvery = 1000 }
        elseif ($dt -eq "1e-5") { $logEvery = 10000 }
        else { $logEvery = 100000 }

        $csv = "taylor_K${K}_dt${dt}_g${gamma}_N${N}_x${xmax}.csv"
        & $exe `
        --gamma $gamma --N $N --xmax $xmax --Umax $Umax `
        --evolve taylor --K $K --dt $dt --tmax $tmax `
        --log $logEvery --csv $csv `
        --no_wide no-wide
    }
}
