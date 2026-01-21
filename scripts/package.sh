#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
app_dir="$root_dir/app"
src_dir="$root_dir/src"
dist_dir="$root_dir/dist"

# Auto-detect platform
detect_platform() {
  local os arch
  os="$(uname -s | tr '[:upper:]' '[:lower:]')"
  arch="$(uname -m)"

  # Normalize OS name
  case "$os" in
    darwin) os="darwin" ;;
    linux)  os="linux" ;;
    *)      echo "Unsupported OS: $os" >&2; exit 1 ;;
  esac

  # Normalize architecture
  case "$arch" in
    x86_64|amd64) arch="x86_64" ;;
    arm64|aarch64) arch="arm64" ;;
    *)      echo "Unsupported architecture: $arch" >&2; exit 1 ;;
  esac

  echo "${os}-${arch}"
}

mkdir -p "$dist_dir"

package_one() {
  local platform="$1"
  local bin_path="$2"
  local staging="$dist_dir/$platform/staging"
  local zip_name="$dist_dir/fmendel-spc-${platform}.zip"

  if [[ ! -f "$bin_path" ]]; then
    echo "Binary not found: $bin_path" >&2
    echo "Run 'make' first to build the binary." >&2
    exit 1
  fi

  rm -rf "$staging"
  mkdir -p "$staging"

  cp "$app_dir"/* "$staging/"
  cp "$bin_path" "$staging/mendel"
  chmod +x "$staging/mendel"

  (cd "$staging" && zip -q -r "$zip_name" .)
  rm -rf "$staging"
  echo "wrote $zip_name"
}

platform="$(detect_platform)"
echo "Detected platform: $platform"
package_one "$platform" "$src_dir/mendel"
