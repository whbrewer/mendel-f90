#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
app_dir="$root_dir/app"
bin_dir="$root_dir/bin"
dist_dir="$root_dir/dist"

mkdir -p "$dist_dir"

package_one() {
  local platform="$1"
  local bin_path="$2"
  local out_dir="$dist_dir/$platform"
  local staging="$out_dir/staging"
  local zip_name="$dist_dir/fmendel-spc-${platform}.zip"

  rm -rf "$staging"
  mkdir -p "$staging"

  cp "$app_dir"/* "$staging/"
  cp "$bin_path" "$staging/mendel"
  chmod +x "$staging/mendel"

  (cd "$staging" && zip -q -r "$zip_name" .)
  rm -rf "$staging"
  echo "wrote $zip_name"
}

package_one "linux-x86_64" "$bin_dir/linux-x86_64/mendel"
package_one "darwin-x86_64" "$bin_dir/darwin-x86_64/mendel"
package_one "darwin-arm64" "$bin_dir/darwin-arm64/mendel"
