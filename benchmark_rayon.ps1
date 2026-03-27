# PowerShell script to benchmark Rayon parallelization performance

Write-Host "Benchmarking Rayon Parallelization Performance"
Write-Host "=============================================="

# Configuration
$iterations = 3
$hpoTerm = "HP:0001083"
$populationSize = 200
$generations = 20

# Function to run GA and measure time
function Measure-GATime {
    param(
        [string]$threadCount,
        [int]$runNumber
    )

    Write-Host "Run $runNumber with $threadCount threads: " -NoNewline

    # Set environment variable
    if ($threadCount -eq "serial") {
        $env:RAYON_NUM_THREADS = "1"
        $label = "Serial (1 thread)"
    } elseif ($threadCount -eq "parallel") {
        $env:RAYON_NUM_THREADS = $null
        $label = "Parallel (auto-detect)"
    } else {
        $env:RAYON_NUM_THREADS = $threadCount
        $label = "$threadCount threads"
    }

    # Run the GA and measure time
    $startTime = Get-Date
    $result = & cargo run --release --bin ga -- --hpo-term $hpoTerm -p $populationSize -g $generations 2>$null
    $endTime = Get-Date

    $duration = ($endTime - $startTime).TotalSeconds

    Write-Host ("{0:F2} seconds" -f $duration)

    # Extract thread count from output
    $threadLine = $result | Where-Object { $_ -match "Rayon threads: (\d+)" }
    if ($threadLine -match "Rayon threads: (\d+)") {
        $actualThreads = $matches[1]
        Write-Host "  Actual threads used: $actualThreads"
    }

    return $duration
}

# Serial runs
Write-Host "`n=== SERIAL RUNS (1 thread) ==="
$serialTimes = @()
for ($i = 1; $i -le $iterations; $i++) {
    $time = Measure-GATime -threadCount "serial" -runNumber $i
    $serialTimes += $time
}

# Parallel runs
Write-Host "`n=== PARALLEL RUNS (auto-detect threads) ==="
$parallelTimes = @()
for ($i = 1; $i -le $iterations; $i++) {
    $time = Measure-GATime -threadCount "parallel" -runNumber $i
    $parallelTimes += $time
}

# Calculate statistics
$serialAvg = ($serialTimes | Measure-Object -Average).Average
$parallelAvg = ($parallelTimes | Measure-Object -Average).Average
$speedup = $serialAvg / $parallelAvg

Write-Host "`n=== RESULTS ==="
Write-Host ("Serial average:   {0:F2} seconds" -f $serialAvg)
Write-Host ("Parallel average: {0:F2} seconds" -f $parallelAvg)
Write-Host ("Speedup:          {0:F2}x" -f $speedup)
Write-Host ("Performance gain: {0:F1}%" -f (($speedup - 1) * 100))

# Clean up
$env:RAYON_NUM_THREADS = $null
