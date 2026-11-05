# Constructive upper bounds on F(n) by local search on edge signs.

function New-Rng {
    param([int]$Seed)
    return [System.Random]::new($Seed)
}

function New-SeidelMatrix {
    param([int]$N, [System.Random]$Rng)
    $a = [int[,]]::new($N, $N)
    for ($i = 0; $i -lt $N; $i++) {
        for ($j = $i + 1; $j -lt $N; $j++) {
            $s = if ($Rng.Next(2) -eq 0) { -1 } else { 1 }
            $a[$i, $j] = $s
            $a[$j, $i] = $s
        }
    }
    return , $a
}

function Get-Score {
    param([int[,]]$A, [int[]]$X)
    $n = $X.Length
    $total = 0
    for ($i = 0; $i -lt $n; $i++) {
        for ($j = $i + 1; $j -lt $n; $j++) {
            $total += $A[$i, $j] * $X[$i] * $X[$j]
        }
    }
    return $total
}

function Get-MaxAbsScoreApprox {
    param(
        [int[,]]$A,
        [System.Random]$Rng,
        [int]$NRestarts = 16,
        [int]$MaxPasses = 40
    )
    $n = $A.GetLength(0)
    $best = 0
    for ($r = 0; $r -lt $NRestarts; $r++) {
        $x = @(0) * $n
        for ($i = 0; $i -lt $n; $i++) {
            $x[$i] = if ($Rng.Next(2) -eq 0) { -1 } else { 1 }
        }
        if ($x[$n - 1] -eq -1) {
            for ($i = 0; $i -lt $n; $i++) { $x[$i] = -$x[$i] }
        }
        $cur = [math]::Abs((Get-Score -A $A -X $x))
        for ($pass = 0; $pass -lt $MaxPasses; $pass++) {
            $improved = $false
            $order = 0..($n - 2) | Sort-Object { $Rng.Next() }
            foreach ($i in $order) {
                $x[$i] = -$x[$i]
                $val = [math]::Abs((Get-Score -A $A -X $x))
                if ($val -gt $cur) {
                    $cur = $val
                    $improved = $true
                } else {
                    $x[$i] = -$x[$i]
                }
            }
            if (-not $improved) { break }
        }
        if ($cur -gt $best) { $best = $cur }
    }
    return [int]$best
}

function Invoke-LocalSearchUpper {
    param(
        [int]$N,
        [int]$Seed = 0,
        [int]$NOuter = 80,
        [int]$NRestarts = 12,
        [double]$Temperature0 = 2.0
    )
    $rng = New-Rng -Seed $Seed
    $a = New-SeidelMatrix -N $N -Rng $rng
    $cur = Get-MaxAbsScoreApprox -A $a -Rng $rng -NRestarts $NRestarts
    $best = $cur
    $bestA = $a.Clone()
    $evals = 1
    $pairs = @()
    for ($i = 0; $i -lt $N; $i++) {
        for ($j = $i + 1; $j -lt $N; $j++) {
            $pairs += , @($i, $j)
        }
    }

    for ($t = 1; $t -le $NOuter; $t++) {
        $temperature = $Temperature0 * (1.0 - $t / ($NOuter + 1))
        $pj = $pairs[$rng.Next($pairs.Count)]
        $i = $pj[0]; $j = $pj[1]
        $a[$i, $j] *= -1
        $a[$j, $i] *= -1
        $val = Get-MaxAbsScoreApprox -A $a -Rng $rng -NRestarts ([math]::Max(4, [int]($NRestarts / 2)))
        $evals++
        $delta = $val - $cur
        $accept = $false
        if ($delta -le 0) {
            $accept = $true
        } else {
            $prob = [math]::Exp(-$delta / [math]::Max($temperature, 1e-9))
            if ($rng.NextDouble() -lt $prob) { $accept = $true }
        }
        if ($accept) {
            $cur = $val
            if ($val -lt $best) {
                $best = $val
                $bestA = $a.Clone()
            }
        } else {
            $a[$i, $j] *= -1
            $a[$j, $i] *= -1
        }
    }

    $final = Get-MaxAbsScoreApprox -A $bestA -Rng $rng -NRestarts ([math]::Max(48, 4 * $NRestarts)) -MaxPasses 80
    $evals++
    $upper = [math]::Max($best, $final)
    return [ordered]@{
        n          = $N
        upper_f    = [int]$upper
        normalized = $upper / [math]::Pow($N, 1.5)
        method     = "sa_edge_flips+greedy_x"
        seed       = $Seed
        evals      = $evals
    }
}

function Test-IsPrime {
    param([int]$N)
    if ($N -lt 2) { return $false }
    if (($N % 2) -eq 0) { return ($N -eq 2) }
    $r = [int][math]::Floor([math]::Sqrt($N))
    for ($d = 3; $d -le $r; $d += 2) {
        if (($N % $d) -eq 0) { return $false }
    }
    return $true
}

function Get-PaleyUpper {
    param([int]$N, [int]$Seed = 0)
    if (($N % 4) -ne 3) { return $null }
    if (-not (Test-IsPrime -N $N)) { return $null }
    $residues = @{}
    for ($k = 1; $k -lt $N; $k++) {
        $residues[(($k * $k) % $N)] = $true
    }
    $a = [int[,]]::new($N, $N)
    for ($i = 0; $i -lt $N; $i++) {
        for ($j = $i + 1; $j -lt $N; $j++) {
            $d = ($j - $i) % $N
            $s = if ($residues.ContainsKey($d)) { 1 } else { -1 }
            $a[$i, $j] = $s
            $a[$j, $i] = $s
        }
    }
    $rng = New-Rng -Seed $Seed
    $upper = Get-MaxAbsScoreApprox -A $a -Rng $rng -NRestarts 96 -MaxPasses 100
    return [ordered]@{
        n          = $N
        upper_f    = [int]$upper
        normalized = $upper / [math]::Pow($N, 1.5)
        method     = "paley_seidel"
        seed       = $Seed
        evals      = 1
    }
}
