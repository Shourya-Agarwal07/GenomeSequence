Write-Host "Checking GenomeToolkit development tools..."
foreach ($tool in @("cl", "cmake", "git", "code")) {
    $found = Get-Command $tool -ErrorAction SilentlyContinue
    if ($found) {
        Write-Host "[FOUND] $tool : $($found.Source)"
    } else {
        Write-Host "[MISSING] $tool"
    }
}
