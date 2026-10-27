#!/usr/bin/env bash
# Build the library versions named in the manuscript into a project-local prefix.
set -euo pipefail
task_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
task_deps="$task_root/.deps"
task_prefix="$task_deps/paper"
task_jobs="${JOBS:-1}"
mkdir -p "$task_deps/downloads" "$task_deps/src" "$task_prefix"

fetch_checked() {
    local filename="$1" url="$2" sha="$3"
    local archive="$task_deps/downloads/$filename"
    if [[ ! -f "$archive" ]]; then
        curl --fail --location --retry 3 "$url" --output "$archive.partial"
        mv "$archive.partial" "$archive"
    fi
    printf '%s  %s\n' "$sha" "$archive" | sha256sum --check --status
}

fetch_checked flint-2.8.0.tar.gz \
    https://flintlib.org/download/flint-2.8.0.tar.gz \
    584235cdc39d779d9920eaef16fe084f3c26ffeeea003a3fff64a20a0f33449e
fetch_checked libsodium-1.0.18.tar.gz \
    https://download.libsodium.org/libsodium/releases/old/libsodium-1.0.18.tar.gz \
    6f504490b342a4f8a4c4a02fc9b866cbef8622d5df4e5452b46be121e46636c1

tar -xzf "$task_deps/downloads/flint-2.8.0.tar.gz" -C "$task_deps/src"
tar -xzf "$task_deps/downloads/libsodium-1.0.18.tar.gz" -C "$task_deps/src"
(
    cd "$task_deps/src/flint-2.8.0"
    ./configure --prefix="$task_prefix" --with-gmp=/usr --with-mpfr=/usr --disable-static
    make -j "$task_jobs"
    make install
)
(
    cd "$task_deps/src/libsodium-1.0.18"
    ./configure --prefix="$task_prefix" --disable-static
    make -j "$task_jobs"
    make check
    make install
)
printf 'export PV_PREFIX=%q\n' "$task_prefix" > "$task_prefix/env.sh"
printf 'export LD_LIBRARY_PATH=%q:${LD_LIBRARY_PATH:-}\n' "$task_prefix/lib" >> "$task_prefix/env.sh"
printf 'export PKG_CONFIG_PATH=%q:${PKG_CONFIG_PATH:-}\n' "$task_prefix/lib/pkgconfig" >> "$task_prefix/env.sh"
printf 'Installed FLINT 2.8.0 and libsodium 1.0.18. Run: source %q\n' "$task_prefix/env.sh"
