# Experiment 2: existing majorants vs vanishing one-sided modulus.
[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot "Exact.ps1")

$ResultsDir = Join-Path $Root "results"
New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null
$t0 = Get-Date

function Get-Alpha([int]$N, [int]$F) {
    return $F / [math]::Pow($N, 1.5)
}

function Get-PaddingAlpha([int]$N, [int]$F, [int]$D) {
    # F(n+d) <= F(n) + d(n-2) + d(d-1)/2   (n>=5; for n<5 still a valid crude upper via n-step)
    $stepBase = [math]::Max($N - 2, 0)
    $upperF = $F + $D * $stepBase + ($D * ($D - 1)) / 2.0
    return $upperF / [math]::Pow($N + $D, 1.5)
}

$ExactF = $script:A399208
Write-Host "== Experiment 2: majorants and neutral increments =="

$series = @()
for ($n = 1; $n -le 15; $n++) {
    $f = $ExactF[$n - 1]
    $a = Get-Alpha $n $f
    $delta = if ($n -lt 15) { $ExactF[$n] - $f } else { $null }
    $hard = if ($n -ge 5 -and $n -lt 15) { $n - 2 } else { $null }
    $neutral = if ($n -lt 15) {
        $f * ([math]::Pow(1.0 + 1.0 / $n, 1.5) - 1.0)
    } else { $null }
    $rho = 0.5 - $a
    $super = $false
    if ($null -ne $delta -and $null -ne $neutral) {
        $super = ([double]$delta -gt [double]$neutral + 1e-12)
    }

    $pMax = $null
    $pArg = $null
    if ($n -ge 5) {
        $best = [double]::NegativeInfinity
        $bestD = 0
        for ($d = 1; $d -le 500; $d++) {
            $p = (Get-PaddingAlpha $n $f $d) - $a
            if ($p -gt $best) { $best = $p; $bestD = $d }
        }
        $pMax = $best
        $pArg = $bestD
    }

    $row = [ordered]@{
        n              = $n
        F              = $f
        alpha_n        = $a
        delta_n        = $delta
        hard_cap_n_minus_2 = $hard
        neutral_r_n    = $neutral
        super_neutral  = $super
        conference_residual = $rho
        padding_excess_max_d500 = $pMax
        padding_excess_argmax_d = $pArg
    }
    $series += $row
    $dn = if ($null -eq $delta) { "   -" } else { "{0,4}" -f $delta }
    $rn = if ($null -eq $neutral) { "   -" } else { "{0:N3}" -f $neutral }
    $pm = if ($null -eq $pMax) { "   -" } else { "{0:N3}" -f $pMax }
    Write-Host ("  n={0,2}  a={1:N6}  dF={2}  r_n={3}  super={4}  rho={5:N6}  Pmax={6}" -f `
        $n, $a, $dn, $rn, $super, $rho, $pm)
}

Write-Host ""
Write-Host "== Hypothesis H2 checks =="

# (1) padding excess >= 0.1 for every n=5..15
$padOk = $true
$padRows = @()
foreach ($r in $series) {
    if ([int]$r.n -lt 5) { continue }
    $ok = [double]$r.padding_excess_max_d500 -ge 0.1 - 1e-15
    $padRows += [ordered]@{ n = $r.n; Pmax = $r.padding_excess_max_d500; argmax_d = $r.padding_excess_argmax_d; ok = $ok }
    if (-not $ok) { $padOk = $false }
}

# (2) conference residual on n=8..15 stays >= 0.03 and not monotone to 0
$residMin = [double]::PositiveInfinity
$residMonoDown = $true
$prev = $null
$residRows = @()
foreach ($r in $series) {
    $n = [int]$r.n
    if ($n -lt 8) { continue }
    $rho = [double]$r.conference_residual
    if ($rho -lt $residMin) { $residMin = $rho }
    if ($null -ne $prev -and $rho -gt $prev + 1e-15) { $residMonoDown = $false }
    $prev = $rho
    $residRows += [ordered]@{ n = $n; rho = $rho }
}
$residFloorOk = ($residMin -ge 0.03 - 1e-15)
$pred2Ok = $residFloorOk -and (-not $residMonoDown)

# (3) at least three super-neutral steps in n=5..14
$superNs = @()
foreach ($r in $series) {
    $n = [int]$r.n
    if ($n -ge 5 -and $n -le 14 -and $r.super_neutral) { $superNs += $n }
}
$pred3Ok = ($superNs.Count -ge 3)

$h2 = [ordered]@{
    prediction_1_padding_excess_ge_0_1 = [ordered]@{
        holds = $padOk
        rows  = $padRows
        note  = "Iterated F(n+1)<=F(n)+n-2 majorant; P_n(d)->infty as d->infty."
    }
    prediction_2_conference_residual_stubborn = [ordered]@{
        holds            = $pred2Ok
        min_rho_n_8_15   = $residMin
        monotone_decreasing = $residMonoDown
        rows             = $residRows
        note             = "epsilon_n = 1/2 - alpha_n does not visibly tend to 0 on the exact prefix."
    }
    prediction_3_at_least_three_super_neutral_steps = [ordered]@{
        holds = $pred3Ok
        count = $superNs.Count
        ns    = $superNs
        note  = "delta_n = F(n+1)-F(n) exceeds constant-alpha neutral increment r_n."
    }
    h2_overall_holds = [bool]($padOk -and $pred2Ok -and $pred3Ok)
}

Write-Host ("  prediction_1_padding_excess_ge_0_1: holds={0}" -f $h2.prediction_1_padding_excess_ge_0_1.holds)
Write-Host ("  prediction_2_conference_residual_stubborn: holds={0}  (min rho={1:N6}, mono_down={2})" -f `
    $h2.prediction_2_conference_residual_stubborn.holds, $residMin, $residMonoDown)
Write-Host ("  prediction_3_at_least_three_super_neutral_steps: holds={0}  count={1}  ns=[{2}]" -f `
    $h2.prediction_3_at_least_three_super_neutral_steps.holds, $superNs.Count, ($superNs -join ","))
Write-Host ("  H2 overall: holds={0}" -f $h2.h2_overall_holds)

$elapsed = ((Get-Date) - $t0).TotalSeconds
$payload = [ordered]@{
    hypothesis = "Iterated sharp-step padding and the conference ceiling 1/2 do not supply a vanishing one-sided all-future modulus on the exact F(1..15) table; several super-neutral one-step increments reopen local excess."
    open_problem = "existence of lim alpha_n (Ivanisvili, MathOverflow 413935)"
    attack_source = "HANDOFF 2026-11-04 + CORE.md one-step bound n-2"
    prior_experiment = "experiment1 (finite E_n under n^{-1/2})"
    series = $series
    h2 = $h2
    implication = "A proof via the one-sided criterion needs a cross-order construction strictly below the conference gap, or a rate alpha_n -> L; recursive n-2 stepping alone cannot close it."
    elapsed_sec = [math]::Round($elapsed, 3)
    limitations = @(
        "Analytic ceiling 1/2 is a limsup bound, not a proved per-n bound.",
        "Exact prefix only through n=15.",
        "Does not prove or disprove convergence."
    )
}

$jsonPath = Join-Path $ResultsDir "experiment2.json"
($payload | ConvertTo-Json -Depth 8) | Set-Content -Path $jsonPath -Encoding UTF8

$csvPath = Join-Path $ResultsDir "experiment2_series.csv"
$csv = New-Object System.Text.StringBuilder
[void]$csv.AppendLine("n,F,alpha_n,delta_n,hard_cap,neutral_r_n,super_neutral,conference_residual,padding_excess_max_d500,padding_argmax_d")
foreach ($r in $series) {
    [void]$csv.AppendLine(("{0},{1},{2:R},{3},{4},{5},{6},{7:R},{8},{9}" -f `
        $r.n, $r.F, $r.alpha_n,
        $(if ($null -eq $r.delta_n) { "" } else { $r.delta_n }),
        $(if ($null -eq $r.hard_cap_n_minus_2) { "" } else { $r.hard_cap_n_minus_2 }),
        $(if ($null -eq $r.neutral_r_n) { "" } else { [string]$r.neutral_r_n }),
        $r.super_neutral,
        $r.conference_residual,
        $(if ($null -eq $r.padding_excess_max_d500) { "" } else { [string]$r.padding_excess_max_d500 }),
        $(if ($null -eq $r.padding_excess_argmax_d) { "" } else { $r.padding_excess_argmax_d })))
}
Set-Content -Path $csvPath -Value $csv.ToString() -Encoding UTF8

Write-Host ""
Write-Host "Wrote $jsonPath"
Write-Host "Wrote $csvPath"
Write-Host ("Elapsed {0}s" -f $payload.elapsed_sec)
