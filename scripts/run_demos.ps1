param(
    [string]$Executable = ""
)

$ErrorActionPreference = "Stop"
$repositoryRoot = Split-Path -Parent $PSScriptRoot

if ([string]::IsNullOrWhiteSpace($Executable)) {
    $releaseExecutable = Join-Path $repositoryRoot "build\Release\genome_tool.exe"
    $defaultExecutable = Join-Path $repositoryRoot "build\genome_tool.exe"
    if (Test-Path -LiteralPath $releaseExecutable) {
        $Executable = $releaseExecutable
    } else {
        $Executable = $defaultExecutable
    }
}

if (-not (Test-Path -LiteralPath $Executable -PathType Leaf)) {
    throw "Executable not found: $Executable. Build the project first."
}

$resolvedExecutable = (Resolve-Path -LiteralPath $Executable).Path
Push-Location $repositoryRoot
try {
    foreach ($scenario in @("normal", "shock", "stress")) {
        Write-Host "`n=== $scenario demonstration ==="
        & $resolvedExecutable --demo $scenario
        if ($LASTEXITCODE -ne 0) {
            throw "The $scenario demonstration failed with exit code $LASTEXITCODE."
        }
    }
} finally {
    Pop-Location
}
