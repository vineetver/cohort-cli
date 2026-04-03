$repo = "vineetver/favor-cli"
$installDir = if ($env:FAVOR_INSTALL_DIR) { $env:FAVOR_INSTALL_DIR } else { "$env:USERPROFILE\.local\bin" }

$latest = (Invoke-RestMethod "https://api.github.com/repos/$repo/releases/latest").tag_name
$url = "https://github.com/$repo/releases/download/$latest/favor-x86_64-windows.zip"

Write-Host "Installing favor $latest..."
Write-Host "  From: $url"
Write-Host "  To:   $installDir\favor.exe"

$tmp = New-TemporaryFile | Rename-Item -NewName { $_.Name + ".zip" } -PassThru
Invoke-WebRequest -Uri $url -OutFile $tmp
New-Item -ItemType Directory -Force -Path $installDir | Out-Null
Expand-Archive -Path $tmp -DestinationPath $installDir -Force
Remove-Item $tmp

if ($env:PATH -notlike "*$installDir*") {
    Write-Host ""
    Write-Host "Add to your PATH:"
    Write-Host "  `$env:PATH += `";$installDir`""
    Write-Host ""
    Write-Host "Or permanently:"
    Write-Host "  [Environment]::SetEnvironmentVariable('PATH', `$env:PATH + ';$installDir', 'User')"
}

Write-Host ""
Write-Host "Installed favor $latest"
Write-Host "Run 'favor setup' to configure."
