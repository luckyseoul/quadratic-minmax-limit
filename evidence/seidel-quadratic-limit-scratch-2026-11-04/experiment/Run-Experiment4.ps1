# Experiment 4: F(16) branches vs one-sided excess.
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot "Exact.ps1")
$ResultsDir = Join-Path $Root "results"
New-Item -ItemType Directory -Force -Path $ResultsDir | Out-Null

function Get-Alpha([int]$N, [int]$F) { $F / [math]::Pow($N, 1.5) }

function Get-SuffixMax([hashtable]$A, [int]$Nmax) {
    $suffix = @{}
    $run = [double]::NegativeInfinity
    for ($n = $Nmax; $n -ge 1; $n--) {
        if ($A[$n] -gt $run) { $run = $A[$n] }
        $suffix[$n] = $run
    }
    return $suffix
}

$base = @{}
for ($n = 1; $n -le 15; $n++) { $base[$n] = Get-Alpha $n $script:A399208[$n - 1] }
$suf15 = Get-SuffixMax $base 15
$e8_15 = $suf15[8] - $base[8]

$branches = @()
foreach ($f16 in @(28, 30)) {
    $A = @{} + $base
    $A[16] = Get-Alpha 16 $f16
    $suf = Get-SuffixMax $A 16
    $e8 = $suf[8] - $A[8]
    $cov = $true
    $rows = @()
    for ($n = 1; $n -le 16; $n++) {
        $e = $suf[$n] - $A[$n]
        $eps = 1.0 / [math]::Sqrt($n)
        if ($n -ge 8 -and $e -gt $eps + 1e-15) { $cov = $false }
        $rows += [ordered]@{ n = $n; alpha = $A[$n]; E = $e; eps = $eps }
    }
    $neutral = [math]::Pow(2.0, 1.5) * 10.0
    $branches += [ordered]@{
        F16 = $f16
        alpha16 = $A[16]
        E8_16 = $e8
        E8_15 = $e8_15
        delta_E8 = $e8 - $e8_15
        worsened = ($e8 -gt $e8_15 + 1e-15)
        covered = $cov
        sub_neutral = ($f16 + 0.0 -lt $neutral)
        series = $rows
    }
    Write-Host ("F16={0} a16={1:N6} E8={2:N6} dE8={3:N6} covered={4} sub={5}" -f `
        $f16, $A[16], $e8, ($e8 - $e8_15), $cov, ($f16 -lt $neutral))
}

$b28 = $branches[0]
$b30 = $branches[1]
$p1 = [bool]((-not $b28.worsened) -and $b28.covered -and $b28.sub_neutral)
$p2 = [bool]($b30.worsened -and ($b30.delta_E8 -gt 0.002) -and $b30.covered -and (-not $b30.sub_neutral))
$p3 = [bool]($b28.sub_neutral -and (-not $b30.sub_neutral))
$overall = [bool]($p1 -and $p2 -and $p3)
Write-Host ("H4 p1={0} p2={1} p3={2} overall={3}" -f $p1, $p2, $p3, $overall)

$payload = [ordered]@{
    hypothesis = "F(16)=28 does not worsen E8 and is sub-neutral; F(16)=30 worsens E8 by >0.002, remains under n^{-1/2}, super-neutral."
    e8_15 = $e8_15
    branches = $branches
    h4 = [ordered]@{ p1 = $p1; p2 = $p2; p3 = $p3; overall = $overall }
    implication = "Finite one-sided tests prefer 28 but do not rule out 30. Next bottleneck: certify F(16)."
}
($payload | ConvertTo-Json -Depth 8) | Set-Content (Join-Path $ResultsDir "experiment4.json") -Encoding UTF8
Write-Host "Wrote results/experiment4.json"
