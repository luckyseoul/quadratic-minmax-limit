# Experiment 3: exact doubling vs constant-alpha neutral composition.
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

$ExactF = $script:A399208
$sqrt2_3 = [math]::Pow(2.0, 1.5)  # 2^{3/2}

Write-Host "== Experiment 3: doubling neutrality =="
Write-Host ("  2^(3/2) = {0:R}" -f $sqrt2_3)

$doublings = @()
foreach ($n in 1..7) {
    $m = 2 * $n
    if ($m -gt 15) { continue }
    $fn = $ExactF[$n - 1]
    $fm = $ExactF[$m - 1]
    $neutral = $sqrt2_3 * $fn
    $an = Get-Alpha $n $fn
    $am = Get-Alpha $m $fm
    $eta = $am - $an
    $sub = $fm + 0.0 -lt $neutral
    $row = [ordered]@{
        n              = $n
        m              = $m
        F_n            = $fn
        F_m            = $fm
        neutral_2n     = $neutral
        slack_F        = $neutral - $fm   # >0 means sub-neutral
        alpha_n        = $an
        alpha_m        = $am
        eta            = $eta
        sub_neutral    = $sub
    }
    $doublings += $row
    Write-Host ("  {0,2}->{1,2}  F={2,2}->{3,2}  neutral={4:N3}  slack={5:N3}  eta={6:N6}  sub={7}" -f `
        $n, $m, $fn, $fm, $neutral, ($neutral - $fm), $eta, $sub)
}

Write-Host ""
Write-Host "== Near-doublings (m/n in [1.8,2.2], n>=5, m<=15) =="
$near = @()
for ($n = 5; $n -le 14; $n++) {
    for ($m = $n + 1; $m -le 15; $m++) {
        $ratio = $m / [double]$n
        if ($ratio -lt 1.8 -or $ratio -gt 2.2) { continue }
        $fn = $ExactF[$n - 1]
        $fm = $ExactF[$m - 1]
        $neutral = $fn * [math]::Pow($ratio, 1.5)
        $an = Get-Alpha $n $fn
        $am = Get-Alpha $m $fm
        $row = [ordered]@{
            n           = $n
            m           = $m
            ratio       = $ratio
            F_n         = $fn
            F_m         = $fm
            neutral_F   = $neutral
            slack_F     = $neutral - $fm
            eta         = $am - $an
            sub_neutral = ($fm + 0.0 -lt $neutral)
        }
        $near += $row
        Write-Host ("  {0,2}->{1,2}  ratio={2:N3}  slackF={3:N3}  eta={4:N6}  sub={5}" -f `
            $n, $m, $ratio, ($neutral - $fm), ($am - $an), $row.sub_neutral)
    }
}

Write-Host ""
Write-Host "== F(16) knife-edge at n=8 =="
$f8 = $ExactF[7]
$neutral16 = $sqrt2_3 * $f8
$a8 = Get-Alpha 8 $f8
$cases = @()
foreach ($f16 in @(28, 30)) {
    $a16 = Get-Alpha 16 $f16
    $eta = $a16 - $a8
    $cases += [ordered]@{
        F16         = $f16
        alpha_16    = $a16
        eta_8       = $eta
        sub_neutral = ($f16 + 0.0 -lt $neutral16)
        slack_F     = $neutral16 - $f16
    }
    Write-Host ("  F(16)={0}  alpha16={1:N6}  eta8={2:N6}  neutral={3:N3}  slack={4:N3}  sub={5}" -f `
        $f16, $a16, $eta, $neutral16, ($neutral16 - $f16), ($f16 -lt $neutral16))
}

Write-Host ""
Write-Host "== Hypothesis H3 checks =="
$subExact = @($doublings | Where-Object { $_.sub_neutral })
$pred1Ok = ($subExact.Count -ge 1)
$mature = @($doublings | Where-Object { [int]$_.n -ge 5 })
$etas = @($mature | ForEach-Object { [double]$_.eta })
$worstEta = if ($etas.Count -gt 0) { ($etas | Measure-Object -Maximum).Maximum } else { $null }
$pred2Ok = ($null -ne $worstEta) -and ($worstEta -le 0.10 + 1e-15)

$case28 = $cases | Where-Object { $_.F16 -eq 28 } | Select-Object -First 1
$case30 = $cases | Where-Object { $_.F16 -eq 30 } | Select-Object -First 1
$pred3Ok = ($case28.sub_neutral -eq $true) -and ($case30.sub_neutral -eq $false)

$h3 = [ordered]@{
    prediction_1_some_exact_subneutral_doubling = [ordered]@{
        holds = $pred1Ok
        pairs = @($subExact | ForEach-Object { "{0}->{1}" -f $_.n, $_.m })
        rows  = $subExact
    }
    prediction_2_worst_mature_doubling_eta_le_0_10 = [ordered]@{
        holds     = $pred2Ok
        from_n    = 5
        worst_eta = $worstEta
        rows      = $mature
    }
    prediction_3_F16_on_neutral_knife_edge = [ordered]@{
        holds              = $pred3Ok
        neutral_F16        = $neutral16
        F8                 = $f8
        case_28_subneutral = $case28.sub_neutral
        case_30_subneutral = $case30.sub_neutral
        cases              = $cases
    }
    h3_overall_holds = [bool]($pred1Ok -and $pred2Ok -and $pred3Ok)
}

Write-Host ("  prediction_1_some_exact_subneutral_doubling: holds={0}  count={1}" -f $pred1Ok, $subExact.Count)
Write-Host ("  prediction_2_worst_mature_doubling_eta_le_0_10: holds={0}  worst_eta={1:N6}" -f $pred2Ok, $worstEta)
Write-Host ("  prediction_3_F16_on_neutral_knife_edge: holds={0}" -f $pred3Ok)
Write-Host ("  H3 overall: holds={0}" -f $h3.h3_overall_holds)

$elapsed = ((Get-Date) - $t0).TotalSeconds
$payload = [ordered]@{
    hypothesis = "Exact doublings already achieve sub-neutral composition; worst mature (n>=5) doubling eta <= 0.10; F(16) in {28,30} decides whether n=8 doubling stays sub-neutral."
    open_problem = "existence of lim alpha_n; also F(16) in {28,30}"
    attack_source = "cross-order / doubling form of HANDOFF 2026-11-04 one-sided criterion"
    doublings = $doublings
    near_doublings = $near
    f16_knife_edge = [ordered]@{
        F8          = $f8
        alpha_8     = $a8
        neutral_F16 = $neutral16
        cases       = $cases
    }
    h3 = $h3
    implication = "Deciding F(16) is the immediate computational bottleneck for the first unavailable doubling test of the one-sided criterion. Exact pair 7->14 is already sub-neutral."
    elapsed_sec = [math]::Round($elapsed, 3)
    limitations = @(
        "Only exact through 15; F(16) open.",
        "Doubling is one ratio; all-future needs every m>=n.",
        "Does not prove convergence."
    )
}

$jsonPath = Join-Path $ResultsDir "experiment3.json"
($payload | ConvertTo-Json -Depth 8) | Set-Content -Path $jsonPath -Encoding UTF8

$csvPath = Join-Path $ResultsDir "experiment3_doublings.csv"
$csv = New-Object System.Text.StringBuilder
[void]$csv.AppendLine("n,m,F_n,F_m,neutral_2n,slack_F,alpha_n,alpha_m,eta,sub_neutral")
foreach ($r in $doublings) {
    [void]$csv.AppendLine(("{0},{1},{2},{3},{4:R},{5:R},{6:R},{7:R},{8:R},{9}" -f `
        $r.n, $r.m, $r.F_n, $r.F_m, $r.neutral_2n, $r.slack_F, $r.alpha_n, $r.alpha_m, $r.eta, $r.sub_neutral))
}
Set-Content -Path $csvPath -Value $csv.ToString() -Encoding UTF8

Write-Host ""
Write-Host "Wrote $jsonPath"
Write-Host "Wrote $csvPath"
Write-Host ("Elapsed {0}s" -f $payload.elapsed_sec)
