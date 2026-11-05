# Exact F(n) via Seidel switching reduction (feasible for small n).

$script:A399208 = @(0, 1, 3, 4, 4, 5, 9, 10, 12, 13, 17, 18, 20, 21, 27)

function Get-ExactF {
    param([int]$N)
    if ($N -le 1) { return 0 }
    if ($N -eq 2) { return 1 }

    $m = $N - 1
    $pairs = @()
    for ($i = 0; $i -lt $m; $i++) {
        for ($j = $i + 1; $j -lt $m; $j++) {
            $pairs += , @($i, $j)
        }
    }
    $nPairs = $pairs.Count
    $best = $null

    $states = @()
    for ($s = 0; $s -lt (1 -shl $m); $s++) {
        $x = @(0) * $m
        for ($i = 0; $i -lt $m; $i++) {
            if ((($s -shr $i) -band 1) -eq 1) { $x[$i] = 1 } else { $x[$i] = -1 }
        }
        $states += , $x
    }

    for ($g = 0; $g -lt (1 -shl $nPairs); $g++) {
        $signs = @(0) * $nPairs
        for ($k = 0; $k -lt $nPairs; $k++) {
            if ((($g -shr $k) -band 1) -eq 1) { $signs[$k] = -1 } else { $signs[$k] = 1 }
        }
        $mx = 0
        foreach ($x in $states) {
            $q = 0
            foreach ($v in $x) { $q += $v }
            for ($k = 0; $k -lt $nPairs; $k++) {
                $i = $pairs[$k][0]
                $j = $pairs[$k][1]
                $q += $signs[$k] * $x[$i] * $x[$j]
            }
            $absq = [math]::Abs($q)
            if ($absq -gt $mx) { $mx = $absq }
            if ($null -ne $best -and $mx -ge $best) { break }
        }
        if ($null -eq $best -or $mx -lt $best) { $best = $mx }
    }
    return [int]$best
}

function Test-AgainstOeis {
    param([int[]]$Ns)
    $out = @{}
    foreach ($n in $Ns) {
        $computed = Get-ExactF -N $n
        $published = $script:A399208[$n - 1]
        $out[$n] = [ordered]@{
            computed = $computed
            oeis     = $published
            match    = ($computed -eq $published)
        }
    }
    return $out
}
