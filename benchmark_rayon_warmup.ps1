# PowerShell script to benchmark Rayon parallelization with proper warm-up handling

Write-Host "Benchmarking Rayon Parallelization Performance (With Warm-up Analysis)"
Write-Host "======================================================================="

# Configuration
$iterations = 3
$hpoTerm = "HP:0001083"
$populationSize = 200
$generations = 20

# Function to run GA and measure time
function Measure-GATime {
    param(
        [string]$threadCount,
        [int]$runNumber,
        [switch]$warmup
    )

    $label = if ($warmup) { "Warm-up" } else { "Run $runNumber" }

    Write-Host "$label with $threadCount threads: " -NoNewline

    # Set environment variable
    if ($threadCount -eq "serial") {
        $env:RAYON_NUM_THREADS = "1"
        $threadLabel = "Serial (1 thread)"
    } elseif ($threadCount -eq "parallel") {
        $env:RAYON_NUM_THREADS = $null
        $threadLabel = "Parallel (auto-detect)"
    } else {
        $env:RAYON_NUM_THREADS = $threadCount
        $threadLabel = "$threadCount threads"
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
        if (-not $warmup) {
            Write-Host "  Actual threads used: $actualThreads"
        }
    }

    return $duration
}

# Warm-up runs (not counted in final results)
Write-Host "`n=== WARM-UP PHASE ==="
Write-Host "Performing warm-up runs to stabilize performance..."
Measure-GATime -threadCount "serial" -runNumber 0 -warmup | Out-Null
Measure-GATime -threadCount "parallel" -runNumber 0 -warmup | Out-Null
Write-Host "Warm-up complete.`n"

# Serial runs
Write-Host "=== SERIAL RUNS (1 thread) ==="
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

# Calculate standard deviations
$serialVariance = 0
foreach ($time in $serialTimes) {
    $serialVariance += [math]::Pow(($time - $serialAvg), 2)
}
$serialStdDev = [math]::Sqrt($serialVariance / $serialTimes.Count)

$parallelVariance = 0
foreach ($time in $parallelTimes) {
    $parallelVariance += [math]::Pow(($time - $parallelAvg), 2)
}
$parallelStdDev = [math]::Sqrt($parallelVariance / $parallelTimes.Count)

Write-Host "`n=== RESULTS ==="
Write-Host ("Serial average:   {0:F2} ± {1:F2} seconds" -f $serialAvg, $serialStdDev)
Write-Host ("Parallel average: {0:F2} ± {1:F2} seconds" -f $parallelAvg, $parallelStdDev)
Write-Host ("Speedup:          {0:F2}x" -f $speedup)
Write-Host ("Performance gain: {0:F1}%" -f (($speedup - 1) * 100))

# Individual run analysis
Write-Host "`n=== INDIVIDUAL RUN ANALYSIS ==="
Write-Host "Serial runs:   $(($serialTimes | ForEach-Object { '{0:F2}' -f $_ }) -join ', ')"
Write-Host "Parallel runs: $(($parallelTimes | ForEach-Object { '{0:F2}' -f $_ }) -join ', ')"

# Performance stability analysis
$serialCoeffVar = ($serialStdDev / $serialAvg) * 100
$parallelCoeffVar = ($parallelStdDev / $parallelAvg) * 100

Write-Host ("`nSerial variability:   {0:F1}% (coefficient of variation)" -f $serialCoeffVar)
Write-Host ("Parallel variability: {0:F1}% (coefficient of variation)" -f $parallelCoeffVar)

# Clean up
$env:RAYON_NUM_THREADS = $null

Write-Host "`n=== CONCLUSION ==="
if ($speedup -gt 1) {
    Write-Host "Parallel execution is FASTER by $([math]::Round(($speedup - 1) * 100, 1))%"
} elseif ($speedup -lt 1) {
    Write-Host "Serial execution is FASTER by $([math]::Round((1 - $speedup) * 100, 1))%"
} else {
    Write-Host "Performance is approximately equal"
}
