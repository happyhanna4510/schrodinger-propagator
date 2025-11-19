# =============================
# run_plots_all.ps1
# Generuje wykresy dla wielu U0, dt i gamma
# =============================

$resultsRoot = ".\results"
$plotsRoot   = ".\plots_out"

# --- KONFIGURACJA ---
# Jakie wartości U0, gamma i dt uwzględniamy
$U0List    = @(-1, 1, 5, -5, -30, -100)      # dodawaj / usuwaj jak chcesz
$gammaList = @(10)                # na razie tylko gamma=10
$dtList    = @("1e-3")  # dla U0=-1 możesz użyć tylko "1e-5", to nie szkodzi jeśli folderów nie będzie

#$dtList    = @("1e-5","1e-4","1e-3")  # dla U0=-1 możesz użyć tylko "1e-5", to nie szkodzi jeśli folderów nie będzie
# =============================

function Ensure-Dir($path) {
    if (-not (Test-Path -LiteralPath $path)) {
        New-Item -ItemType Directory -Path $path | Out-Null
    }
}

foreach ($U0 in $U0List) {
    $u0FolderName = "U0_$U0"

    foreach ($gamma in $gammaList) {
        $gFolderName = "g$gamma"

        foreach ($dt in $dtList) {
            $dtFolderName = "dt$dt"

            $rootDir = Join-Path $resultsRoot (Join-Path $u0FolderName (Join-Path $gFolderName $dtFolderName))
            $outDir  = Join-Path $plotsRoot   (Join-Path $u0FolderName (Join-Path $gFolderName $dtFolderName))

            if (-not (Test-Path -LiteralPath $rootDir)) {
                Write-Host "Skip: results not found for U0=$U0, gamma=$gamma, dt=$dt  ($rootDir)" -ForegroundColor DarkYellow
                continue
            }

            Write-Host ""
            Write-Host "=== Plots for U0=$U0, gamma=$gamma, dt=$dt ===" -ForegroundColor Cyan
            Write-Host "  root : $rootDir"
            Write-Host "  out  : $outDir"

            Ensure-Dir $outDir

            # wykresy błędów (norm_err, theta_abs)
            python .\scripts\plot_compare_methods.py --root $rootDir --out $outDir

            # heatmap + przekroje |psi|^2
            python .\scripts\plot_abs2_wide.py      --results-dir $rootDir --out-dir $outDir
        }
    }
}

Write-Host ""
Write-Host "=== DONE (run_plots_all.ps1) ==="
