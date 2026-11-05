# Experiment 1: finite diagnostic for the one-sided all-future criterion.
[CmdletBinding()]
param(
    [int]$MaxExact = 6,
    [switch]$SkipExactVerify
)

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot "Exact.ps1")

$ResultsDir = Join-Path $Root "results"
New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null
$t0 = Get-Date

function Get-Alpha([int]$N, [int]$F) {
    if ($N -le 0) { return 0.0 }
    return $F / [math]::Pow($N, 1.5)
}

# Published exact F(n) = A399208; upper envelopes from research archives.
$ExactF = $script:A399208  # indices 0..14 for n=1..15
$UpperF = @{
    16 = 30   # framed Hadamard / both archives
    17 = 32   # Paley puncture + parity
}

Write-Host "== Exact verification vs OEIS A399208 =="
$verifiedForJson = [ordered]@{}
if (-not $SkipExactVerify) {
    $verifyNs = 1..([math]::Max(1, $MaxExact))
    $verified = Test-AgainstOeis -Ns $verifyNs
    foreach ($n in ($verified.Keys | Sort-Object)) {
        $row = $verified[$n]
        $status = if ($row.match) { "OK" } else { "MISMATCH" }
        Write-Host ("  n={0,2}  computed={1,3}  oeis={2,3}  [{3}]" -f $n, $row.computed, $row.oeis, $status)
        $verifiedForJson["$n"] = $row
    }
    if (@($verified.Values | Where-Object { -not $_.match }).Count -gt 0) {
        throw "Exact verification failed; aborting."
    }
} else {
    Write-Host "  (skipped)"
}

Write-Host ""
Write-Host "== Exact alpha_n and one-sided excess E_n(15) =="
$alpha = @{}
$exactRows = @()
for ($n = 1; $n -le 15; $n++) {
    $f = $ExactF[$n - 1]
    $a = Get-Alpha $n $f
    $alpha[$n] = $a
}

# Precompute suffix maxima of alpha on [k,15]
$suffixMax = @{}
$running = [double]::NegativeInfinity
for ($n = 15; $n -ge 1; $n--) {
    if ($alpha[$n] -gt $running) { $running = $alpha[$n] }
    $suffixMax[$n] = $running
}

$rows = @()
for ($n = 1; $n -le 15; $n++) {
    $eExact = $suffixMax[$n] - $alpha[$n]
    $a16 = 30 / [math]::Pow(16, 1.5)
    $a17 = 32 / [math]::Pow(17, 1.5)
    $eBar = [math]::Max($eExact, [math]::Max($a16 - $alpha[$n], [math]::Max($a17 - $alpha[$n], 0.5 - $alpha[$n])))
    $sched = if ($n -ge 1) { 1.0 / [math]::Sqrt($n) } else { [double]::PositiveInfinity }
    $row = [ordered]@{
        n                = $n
        F                = $ExactF[$n - 1]
        alpha_n          = $alpha[$n]
        suffix_max_alpha = $suffixMax[$n]
        E_n_15           = $eExact
        E_bar_n          = $eBar
        eps_inv_sqrt     = $sched
        E_n_over_eps     = if ($sched -gt 0) { $eExact / $sched } else { $null }
        E_bar_over_eps   = if ($sched -gt 0) { $eBar / $sched } else { $null }
        analytic_lower   = if ($n -le 1) { 0.0 } else { [math]::Sqrt(1.0 - 1.0 / $n) / [math]::PI }
        analytic_upper   = 0.5
    }
    $rows += $row
    Write-Host ("  n={0,2}  F={1,3}  a={2:N6}  E_15={3:N6}  Ebar={4:N6}  1/sqrt(n)={5:N6}" -f `
        $n, $row.F, $row.alpha_n, $row.E_n_15, $row.E_bar_n, $row.eps_inv_sqrt)
}

Write-Host ""
Write-Host "== Hypothesis H1 checks =="
$byN = @{}
foreach ($r in $rows) { $byN[[int]$r.n] = $r }

# (1) E_n(15) eventually decreasing for n>=8
$deltas = @()
$monoOk = $true
for ($n = 8; $n -le 14; $n++) {
    $d = [double]$byN[$n + 1].E_n_15 - [double]$byN[$n].E_n_15
    $deltas += [ordered]@{ from = $n; to = ($n + 1); delta_E = $d }
    if ($d -gt 1e-15) { $monoOk = $false }
}

# (2) E_n(15) <= n^{-1/2} for n>=8
$schedExactOk = $true
$schedExactRows = @()
for ($n = 8; $n -le 15; $n++) {
    $ok = [double]$byN[$n].E_n_15 -le [double]$byN[$n].eps_inv_sqrt + 1e-15
    $schedExactRows += [ordered]@{ n = $n; E = $byN[$n].E_n_15; eps = $byN[$n].eps_inv_sqrt; ok = $ok }
    if (-not $ok) { $schedExactOk = $false }
}

# (3) E_bar_n <= n^{-1/2} for n>=8
$schedBarOk = $true
$schedBarRows = @()
for ($n = 8; $n -le 15; $n++) {
    $ok = [double]$byN[$n].E_bar_n -le [double]$byN[$n].eps_inv_sqrt + 1e-15
    $schedBarRows += [ordered]@{ n = $n; Ebar = $byN[$n].E_bar_n; eps = $byN[$n].eps_inv_sqrt; ok = $ok }
    if (-not $ok) { $schedBarOk = $false }
}

$h1 = [ordered]@{
    prediction_1_E_n_eventually_decreasing = [ordered]@{
        from_n = 8
        holds  = $monoOk
        deltas = $deltas
        note   = "Expected fragile: alpha_n oscillates even though m_n is monotone (CORE.md §2)."
    }
    prediction_2_exact_excess_below_inv_sqrt = [ordered]@{
        from_n = 8
        holds  = $schedExactOk
        rows   = $schedExactRows
    }
    prediction_3_mixed_excess_below_inv_sqrt = [ordered]@{
        from_n = 8
        holds  = $schedBarOk
        rows   = $schedBarRows
        note   = "Uses F(16)<=30, F(17)<=32, and analytic ceiling 1/2."
    }
    h1_overall_holds = [bool]($monoOk -and $schedExactOk -and $schedBarOk)
}

Write-Host ("  prediction_1_E_n_eventually_decreasing: holds={0}" -f $h1.prediction_1_E_n_eventually_decreasing.holds)
Write-Host ("  prediction_2_exact_excess_below_inv_sqrt: holds={0}" -f $h1.prediction_2_exact_excess_below_inv_sqrt.holds)
Write-Host ("  prediction_3_mixed_excess_below_inv_sqrt: holds={0}" -f $h1.prediction_3_mixed_excess_below_inv_sqrt.holds)
Write-Host ("  H1 overall: holds={0}" -f $h1.h1_overall_holds)

$elapsed = ((Get-Date) - $t0).TotalSeconds
$payload = [ordered]@{
    hypothesis = "Finite one-sided excess E_n(15) is eventually decreasing for n>=8 and dominated by n^{-1/2}; mixed excess with F(16)<=30, F(17)<=32, alpha<=1/2 likewise."
    open_problem = "existence of lim alpha_n (Ivanisvili, MathOverflow 413935)"
    attack_source = "luckyseoul/quadratic-minmax-limit HANDOFF 2026-11-04 one-sided all-future criterion"
    references = @(
        "https://github.com/luckyseoul/quadratic-minmax-limit",
        "https://github.com/Robby955/mo-413935-research",
        "https://oeis.org/A399208",
        "https://mathoverflow.net/questions/413935"
    )
    exact_verification = $verifiedForJson
    series = $rows
    upper_envelopes = [ordered]@{
        "16" = [ordered]@{ F_upper = 30; alpha_upper = (30 / [math]::Pow(16, 1.5)) }
        "17" = [ordered]@{ F_upper = 32; alpha_upper = (32 / [math]::Pow(17, 1.5)) }
        analytic_ceiling = 0.5
    }
    h1 = $h1
    A399208 = $script:A399208
    elapsed_sec = [math]::Round($elapsed, 3)
    limitations = @(
        "Finite N=15 cannot prove an all-future epsilon_n -> 0.",
        "Mixed excess uses possibly loose upper envelopes.",
        "Does not settle existence of lim alpha_n."
    )
}

$jsonPath = Join-Path $ResultsDir "experiment1.json"
($payload | ConvertTo-Json -Depth 8) | Set-Content -Path $jsonPath -Encoding UTF8

$csvPath = Join-Path $ResultsDir "experiment1_series.csv"
$csv = New-Object System.Text.StringBuilder
[void]$csv.AppendLine("n,F,alpha_n,suffix_max_alpha,E_n_15,E_bar_n,eps_inv_sqrt,E_n_over_eps,E_bar_over_eps,analytic_lower,analytic_upper")
foreach ($r in $rows) {
    [void]$csv.AppendLine(("{0},{1},{2:R},{3:R},{4:R},{5:R},{6:R},{7:R},{8:R},{9:R},{10:R}" -f `
        $r.n, $r.F, $r.alpha_n, $r.suffix_max_alpha, $r.E_n_15, $r.E_bar_n, `
        $r.eps_inv_sqrt, $r.E_n_over_eps, $r.E_bar_over_eps, $r.analytic_lower, $r.analytic_upper))
}
Set-Content -Path $csvPath -Value $csv.ToString() -Encoding UTF8

Write-Host ""
Write-Host "Wrote $jsonPath"
Write-Host "Wrote $csvPath"
Write-Host ("Elapsed {0}s" -f $payload.elapsed_sec)
