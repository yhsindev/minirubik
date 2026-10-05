# ripes-measure.ps1 - run one Ripes CLI simulation and report
# retired instructions, model execution time, and peak host memory.
#
# Usage:
#   .\ripes-measure.ps1 -Src loop.s -Proc RV32_ISS [-Runs 3]
#
# Peak memory is the Ripes process's PeakWorkingSet64, sampled until exit.
param(
    [Parameter(Mandatory)] [string] $Src,
    [string] $Proc = "RV32_ISS",
    [int]    $Runs = 1,
    [string] $Ripes = "C:\Users\itlab\Desktop\1151\Ripes-v2.2.6-106-g5b8a616-win-x86_64\Ripes.exe",
    [int]    $TimeoutMs = 3600000
)

$srcPath = (Resolve-Path $Src).Path
$results = @()

for ($i = 1; $i -le $Runs; $i++) {
    $out = New-TemporaryFile
    $ripesArgs = @("--mode", "cli", "--src", "`"$srcPath`"", "-t", "asm",
              "--proc", $Proc, "--iret", "--exectime", "--json",
              "--timeout", $TimeoutMs,
              "--reginit", "gpr:2=0x7FFFFFF0,3=0x10000000")
    $p = Start-Process -FilePath $Ripes -ArgumentList $ripesArgs -NoNewWindow `
                       -PassThru -RedirectStandardOutput $out.FullName
    $peak = 0
    while (-not $p.HasExited) {
        try { $p.Refresh(); if ($p.PeakWorkingSet64 -gt $peak) { $peak = $p.PeakWorkingSet64 } } catch {}
        Start-Sleep -Milliseconds 20
    }
    $text = Get-Content $out.FullName -Raw
    Remove-Item $out.FullName
    $json = $null
    if ($text -match '(?s)(\{.*\})') { $json = $Matches[1] | ConvertFrom-Json }
    if ($null -eq $json) { Write-Error "Ripes produced no report:`n$text"; exit 1 }

    $iret = [int64] $json.'# instructions retired'
    $ms   = [int64] $json.'execution time (ms)'
    $results += [pscustomobject]@{
        Run       = $i
        Proc      = $Proc
        Iret      = $iret
        ExecMs    = $ms
        IPS       = if ($ms -gt 0) { [math]::Round($iret * 1000.0 / $ms) } else { $null }
        PeakBytes = $peak
        PeakMiB   = [math]::Round($peak / 1MB, 2)
    }
}
$results | Format-Table -AutoSize
